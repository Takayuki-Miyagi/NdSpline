#include "NdSpline.hh"

NdSpline::~NdSpline()
{}

// constructor of spline
// inputs:
// ks : integer vector for the order of polymial for each abscissa
// ns : integer vector array for the dimension of each abscissa
// x : double vector for the points of each abscissa. One can use the concatenated array [x,y,...]
//           x, y, ... are non dcreasing one dimensional array
// f : double vector for the points of each abscissa. One can use fs = reshape( f, (nx*ny*...) )
//         Here f is multidimenstional array f(i,j,...). i,j,... have to be correspond to x,y,..., respectively.
NdSpline::NdSpline(Eigen::VectorXi ks, Eigen::VectorXi ns, Eigen::VectorXd x, Eigen::VectorXd f)
{
  if(ks.size() != ns.size()){
    std::cout << "Error in constructor of spline: Dimension cannot be detected properly!" << std::endl;
    exit(0);
  }

  int n_f = 1; // total number of points
  for(int n=0; n<ns.size(); n++){
    n_f *= ns(n);
  }
  if(f.size() != n_f) {
    std::cout << "Error in constructor of spline: f cannot be detected properly" << std::endl;
    exit(0);
  }

  ndim = ks.size();
  nf_in = n_f;
  int n_start = 0;
  for(int n=0; n<ndim; n++){
    if(ks(n) < 1){
      std::cout << "Error in constructor of spline: polynomial order k has to be greater than 0" << std::endl;
      exit(0);
    }
    Ndrctn.emplace_back(ks(n), Eigen::VectorXd(x.segment(n_start, ns(n))));
    n_start += ns(n);
  }
  coefs.resize(f.size());
  coefs = f;
  set_coefs();
}

// interpolation using the matrix-matrix multiplication
// inputs:
// ns: integer vector specifying the dimenstion
// x:  double vector array specifying the interpolants
void NdSpline::Interpolate(const Eigen::VectorXi& ns, const Eigen::VectorXd& x)
{
  grid_interpolant grid = grid_interpolant(ns, x);
  if(ndim != grid.ndim) {
    std::cout << "Error in interpolate of spline: Dimension would be wrong!" << std::endl;
    return;
  }

  int n = nf_in;
  Eigen::Map<Eigen::MatrixXd>* tmp_org = nullptr;
  Eigen::MatrixXd tmp_int;
  Eigen::MatrixXd coefs_ = coefs;

  for(int i=0; i<ndim; i++){
    if(i==0){
      tmp_org = new Eigen::Map<Eigen::MatrixXd>(coefs_.data(), Ndrctn[i].nx, n / Ndrctn[i].nx);
    }
    else{
      tmp_org = new Eigen::Map<Eigen::MatrixXd>(tmp_int.data(), Ndrctn[i].nx, n / Ndrctn[i].nx);
    }
    interpolate_1d(Ndrctn[i], grid.Ndrctn[i].x, *tmp_org, tmp_int);
    //std::cout << "tmp_org" << std::endl;
    //std::cout << tmp_org->transpose() << std::endl;
    //std::cout << "tmp_int" << std::endl;
    //std::cout << tmp_int << std::endl;
    n = (n * grid.Ndrctn[i].n_interpolant) / Ndrctn[i].nx;
  }
  Eigen::Map<Eigen::VectorXd> tmp_vec(tmp_int.data(), tmp_int.size());
  fs = tmp_vec;
}

void NdSpline::interpolate_1d(const abscissa& abs, const Eigen::VectorXd& x_out, const Eigen::MatrixXd& coefs, Eigen::MatrixXd& f)
{
  int n = abs.nx;
  int m = x_out.size();
  Eigen::MatrixXd BMat(m, n);
#pragma omp parallel for schedule(dynamic)
  for(int j = 0; j < m; j++){
    for(int i = 0; i < n; i++){
      BMat(j,i) = abs.basis_function(i, abs.k, x_out(j));
    }
  }
  f = coefs.transpose() * BMat.transpose();
  //f = BMat * coefs;
  //f = coefs * BMat.transpose();
}

NdSpline::abscissa::abscissa(int k, const Eigen::VectorXd& x)
  : k(k), x(x)
{
  nx = x.size();
  nt = k + nx;
  t.resize(nt);
  set_knot_vector();
//  std::cout << std::fixed << x.transpose() << std::endl;
//  std::cout << std::fixed << t.transpose() << std::endl;
}

void NdSpline::abscissa::set_knot_vector()
{
  t.setConstant(0.0);
  t.segment(0,k).setConstant(x(0));
  t.segment(nx, k).setConstant(x(x.size()-1));

  if(k%2 == 0){
    t.segment(k,nx-k) = x.segment(k/2, nx-k);
    return;
  }

  if(k%2 == 1){
    Eigen::VectorXd tmp(nx-k);
    for(int i=0; i<nx-k; i++){
      tmp(i) = 0.5 * (x(k/2+i) + x(k/2+i+1));
    }
    t.segment(k,nx-k) = tmp;
    return;
  }
}

void NdSpline::set_coefs()
{
  Eigen::MatrixXd tmp, tmp_t;
  for(int i = 0; i<ndim; i++) {
    Eigen::Map<Eigen::MatrixXd>* tmp_coefs = nullptr;
    if(i==0){
      tmp_coefs = new Eigen::Map<Eigen::MatrixXd>(coefs.data(), Ndrctn[i].nx, nf_in / Ndrctn[i].nx);
    }
    else{
      tmp_t = tmp.transpose(); // Eigen is not great
      tmp_coefs = new Eigen::Map<Eigen::MatrixXd>(tmp_t.data(), Ndrctn[i].nx, nf_in / Ndrctn[i].nx);
    }
    tmp.resize(Ndrctn[i].nx, nf_in/Ndrctn[i].nx);
    set_coefs_1d(Ndrctn[i], *tmp_coefs, tmp);
    //std::cout << "tmp_coefs" << std::endl;
    //std::cout << *tmp_coefs << std::endl;
    //std::cout << "tmp" << std::endl;
    //std::cout << tmp << std::endl;
  }
  tmp_t = tmp.transpose();
  Eigen::Map<Eigen::VectorXd> tmp_vec(tmp_t.data(), tmp.size());
  coefs = tmp_vec;
  //std::cout << coefs << std::endl;
}

void NdSpline::set_coefs_1d(const abscissa& abs, const Eigen::MatrixXd& f, Eigen::MatrixXd& coefs)
{
  int nx = abs.nx;
  Eigen::MatrixXd BMat(nx, nx);
#pragma omp parallel for schedule(dynamic)
  for(int j=0; j<nx; j++){
    for(int i=0; i<nx; i++){
      BMat(j,i) = abs.basis_function(i, abs.k, abs.x(j));
    }
  }
  // solve BMat * coefs = f, we want to know coefs
  Eigen::FullPivLU<Eigen::MatrixXd> lu(BMat);
  if (!lu.isInvertible()) {
    std::cerr << "Matrix is singular!" << std::endl;
    exit(0);
  }
  coefs = lu.solve(f);
}

double NdSpline::abscissa::basis_function(int i, int k, double x) const
{
  double f  = this->basis_function_rec(i, k, this->find_interval(x), x);
  return f;
}

//
// basis function phi^{k}_{i}(x) using the Cox-de Boor recursion formula
//
double NdSpline::abscissa::basis_function_rec(int i, int k, int xi, double x) const
{
  if(i+k > t.size()) return 0.0;
  if(k == 1){
    if(i==xi) {
      return 1.0;
    }
    else {
      return 0.0;
    }
  }
  double c1 = (x - t(i)) / (t(i+k-1) - t(i));
  double c2 = (t(i+k)-x) / (t(i+k) - t(i+1));
  if( abs(t(i+k-1) - t(i)) < 1.e-8 ) c1 = 0.0;
  if( abs(t(i+k) - t(i+1)) < 1.e-8 ) c2 = 0.0;
  if( abs(t(i+k-1) - t(i)) < 1.e-8 and abs(x - t(i)) < 1.e-8 ) c1 = 1.0;
  if( abs(t(i+k) - t(i+1)) < 1.e-8 and abs(x - t(i+k)) < 1.e-8 ) c2 = 1.0;
  double f = this->basis_function_rec(i, k-1, xi, x) * c1 + this->basis_function_rec(i+1, k-1, xi, x) * c2;
  return f;
}

//
// return the index i so t(i) <= x <= t(i+1)
// Note that x has to satisfy t(1) <= x <= t(n+k), when the extrapolation is turned off.
//
int NdSpline::abscissa::find_interval(double x) const
{
  if(extrapolation) {
    if(x <= t(k-1)) return k-1;
    if(t(nx-1) <= x) return nx-1;
  }
  else {
    if(x < t.minCoeff() or x > t.maxCoeff()) {
      std::cout << "Error: at line " << __LINE__ << " in " << __FILE__ << std::endl;
      exit(0);
    }
    if(t(k-1) <= x and x <= t(k)) return k-1;
    if(t(nx) <= x and x <= t(nx+1)) return nx;
  }

  for(int i=0; i<nt-1; i++) {
    if(t(i) <= x and x <= t(i+1)) return i;
  }
  std::cout << "Error: at line " << __LINE__ << " in " << __FILE__ << std::endl;
  exit(0);
  return 0;
}

NdSpline::grid_interpolant::grid_interpolant(const Eigen::VectorXi& ns, const Eigen::VectorXd& x)
{
  ndim = ns.size();
  int n_start = 0;
  for(int i=0; i<ndim; i++){
    this->Ndrctn.emplace_back(Eigen::VectorXd(x.segment(n_start, ns(i))));
    n_start += ns(i);
  }
  n_interpolant = 1;
  for(int i=0; i<ndim; i++){
    n_interpolant += Ndrctn[i].n_interpolant;
  }
}

NdSpline::abscissa_interpolant::abscissa_interpolant(const Eigen::VectorXd& x)
  : x(x)
{
  n_interpolant = x.size();
}
