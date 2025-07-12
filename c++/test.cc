#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include "NdSpline.hh"

void test_1d()
{
  double xmin = -2.0, xmax = 2.0;
  int n, k = 4;
  auto f1 = [](double x){ return std::cos(x);};

  Eigen::Vector<double, Eigen::Dynamic> x, f;
  n = 50;
  x.resize(n);
  f.resize(n);
  for (int i=0; i<n; i++){
    x(i) = xmin + (double)(i) / (double)(n-1) * (xmax-xmin);
  }

  for (int i=0; i<n; i++){
    f(i) = f1(x(i));
  }

  Eigen::VectorXi ks(1), ns(1);
  ks << k;
  ns << n;
  NdSpline sp = NdSpline(ks, ns, x, f);

  n = 100;
  x.resize(n);
  f.resize(n);
  for(int i=0; i<n; i++){
    x(i) = xmin + (double)i / (double)(n-1) * (xmax-xmin);
  }

  ns << n;
  sp.Interpolate(ns, x);
  f = sp.fs;

  double max_err = 0.0;
  std::ofstream outfile("file_1d.dat");
  outfile << "x, interpolated, true value" << std::endl;
  for(int i=0; i<n; i++){
    outfile << std::fixed << std::setw(8) << std::setprecision(4) << x(i)
      << std::scientific << std::setw(18) << std::setprecision(6) << f(i)
      << std::scientific << std::setw(18) << std::setprecision(6) << f1(x(i))
      << std::endl;
    max_err = std::max(max_err, (std::abs(sp.fs(i) - f1(x(i)))) / f1(x(i)));
  }
  std::cout << "Maximum interpolation error: " << std::scientific << std::setw(18) << std::setprecision(6) << max_err << std::endl;
}

void test_2d()
{
  double xmin = -2.0, xmax = 2.0;
  double ymin = -1.0, ymax = 3.0;
  int nx=15;
  int ny=16;
  int kx=4;
  int ky=4;
  auto f2 = [](double x, double y){ return std::cos(x)*std::cos(y);};
  //auto f2 = [](double x, double y){ return std::cos(x)*y*y*(y+1);};

  Eigen::Vector<double, Eigen::Dynamic> x, y;
  x.resize(nx);
  y.resize(ny);
  for (int i=0; i<nx; i++){
    x(i) = xmin + (double)(i) / (double)(nx-1) * (xmax-xmin);
  }
  for (int i=0; i<ny; i++){
    y(i) = ymin + (double)(i) / (double)(ny-1) * (ymax-ymin);
  }

  Eigen::Vector<double, Eigen::Dynamic> f;
  f.resize(nx*ny);
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      int idx = j*nx + i;
      f(idx) = f2(x(i), y(j));
    }
  }
  Eigen::VectorXi ks(2), ns(2);
  ks << kx, ky;
  ns << nx, ny;
  Eigen::VectorXd xy(x.size() + y.size());
  xy << x, y;
  NdSpline sp = NdSpline(ks, ns, xy, f);

  nx = 20;
  ny = 20;
  x.resize(nx);
  y.resize(ny);
  f.resize(nx*ny);
  xy.resize(nx+ny);
  for (int i=0; i<nx; i++){
    x(i) = xmin + (double)(i) / (double)(nx-1) * (xmax-xmin);
  }
  for (int i=0; i<ny; i++){
    y(i) = ymin + (double)(i) / (double)(ny-1) * (ymax-ymin);
  }
  ns << nx, ny;
  xy << x, y;
  sp.Interpolate(ns, xy);

  double max_err = 0.0;
  std::ofstream outfile("file_2d.dat");
  outfile << "x, y, interpolated, true value" << std::endl;
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      int idx = j*nx + i;
      outfile << std::fixed << std::setw(8) << std::setprecision(4) << x(i)
        << std::fixed << std::setw(8) << std::setprecision(4) << y(j)
        << std::scientific << std::setw(18) << std::setprecision(6) << sp.fs(idx)
        << std::scientific << std::setw(18) << std::setprecision(6) << f2(x(i), y(j))
        << std::endl;
      max_err = std::max(max_err, (std::abs(sp.fs(idx) - f2(x(i), y(j)))) / f2(x(i), y(j)));
    }
  }
  std::cout << "Maximum interpolation error: " << std::scientific << std::setw(18) << std::setprecision(6) << max_err << std::endl;
}

void test_3d()
{
  double xmin = -2.0, xmax = 2.0;
  double ymin = -1.0, ymax = 3.0;
  double zmin = -3.0, zmax = 1.0;
  int nx=14;
  int ny=15;
  int nz=16;
  int kx=4;
  int ky=4;
  int kz=4;
  //auto f3 = [](double x, double y, double z){ return std::cos(x)*std::cos(y)*z*z;};
  auto f3 = [](double x, double y, double z){ return std::cos(x)*std::cos(y)*cos(z);};

  Eigen::Vector<double, Eigen::Dynamic> x, y, z;
  x.resize(nx);
  y.resize(ny);
  z.resize(nz);
  for (int i=0; i<nx; i++){
    x(i) = xmin + (double)(i) / (double)(nx-1) * (xmax-xmin);
  }
  for (int i=0; i<ny; i++){
    y(i) = ymin + (double)(i) / (double)(ny-1) * (ymax-ymin);
  }
  for (int i=0; i<nz; i++){
    z(i) = zmin + (double)(i) / (double)(nz-1) * (zmax-zmin);
  }

  Eigen::Vector<double, Eigen::Dynamic> f;
  f.resize(nx*ny*nz);
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      for(int k=0; k<nz; k++){
        int idx = k * nx * ny + j * nx + i;
        f(idx) = f3(x(i), y(j), z(k));
      }
    }
  }
  Eigen::VectorXi ks(3), ns(3);
  ks << kx, ky, kz;
  ns << nx, ny, nz;
  Eigen::VectorXd xyz(x.size() + y.size() + z.size());
  xyz << x, y, z;
  NdSpline sp = NdSpline(ks, ns, xyz, f);

  nx = 20;
  ny = 20;
  nz = 20;
  x.resize(nx);
  y.resize(ny);
  z.resize(nz);
  f.resize(nx*ny*nz);
  xyz.resize(nx+ny+nz);
  for (int i=0; i<nx; i++){
    x(i) = xmin + (double)(i) / (double)(nx-1) * (xmax-xmin);
  }
  for (int i=0; i<ny; i++){
    y(i) = ymin + (double)(i) / (double)(ny-1) * (ymax-ymin);
  }
  for (int i=0; i<nz; i++){
    z(i) = zmin + (double)(i) / (double)(nz-1) * (zmax-zmin);
  }
  ns << nx, ny, nz;
  xyz << x, y, z;
  sp.Interpolate(ns, xyz);

  double max_err = 0.0;
  std::ofstream outfile("file_3d.dat");
  outfile << "x, y, z, interpolated, true value" << std::endl;
  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      for(int k=0; k<nz; k++){
        int idx = k * nx * ny + j * nx + i;
        outfile << std::fixed << std::setw(8) << std::setprecision(4) << x(i)
          << std::fixed << std::setw(8) << std::setprecision(4) << y(j)
          << std::fixed << std::setw(8) << std::setprecision(4) << z(k)
          << std::scientific << std::setw(18) << std::setprecision(6) << sp.fs(idx)
          << std::scientific << std::setw(18) << std::setprecision(6) << f3(x(i), y(j), z(k))
          << std::endl;
        max_err = std::max(max_err, (std::abs(sp.fs(idx) - f3(x(i), y(j), z(k)))) / f3(x(i), y(j), z(k)));
      }
    }
  }
  std::cout << "Maximum interpolation error: " << std::scientific << std::setw(18) << std::setprecision(6) << max_err << std::endl;
}

void test_4d()
{
  //double wmin = -5.0, wmax = 3.0;
  //double xmin = -2.0, xmax = 2.0;
  //double ymin = -1.0, ymax = 3.0;
  //double zmin = -3.0, zmax = 1.0;
  double wmin = -2.0, wmax = 2.0;
  double xmin = -2.0, xmax = 2.0;
  double ymin = -2.0, ymax = 2.0;
  double zmin = -2.0, zmax = 2.0;
  int nw=8;
  int nx=8;
  int ny=8;
  int nz=8;
  int kw=4;
  int kx=4;
  int ky=4;
  int kz=4;
  //auto f4 = [](double w, double x, double y, double z){ return std::sin(w)*std::cos(x)*std::cos(y)*z*z;};
  auto f4 = [](double w, double x, double y, double z){ return std::cos(w)*std::cos(x)*std::cos(y)*std::cos(z);};

  Eigen::Vector<double, Eigen::Dynamic> w, x, y, z;
  w.resize(nw);
  x.resize(nx);
  y.resize(ny);
  z.resize(nz);
  for (int i=0; i<nw; i++){
    w(i) = wmin + (double)(i) / (double)(nw-1) * (wmax-wmin);
  }
  for (int i=0; i<nx; i++){
    x(i) = xmin + (double)(i) / (double)(nx-1) * (xmax-xmin);
  }
  for (int i=0; i<ny; i++){
    y(i) = ymin + (double)(i) / (double)(ny-1) * (ymax-ymin);
  }
  for (int i=0; i<nz; i++){
    z(i) = zmin + (double)(i) / (double)(nz-1) * (zmax-zmin);
  }

  Eigen::Vector<double, Eigen::Dynamic> f;
  f.resize(nw*nx*ny*nz);
  for(int h=0; h<nw; h++){
    for(int i=0; i<nx; i++){
      for(int j=0; j<ny; j++){
        for(int k=0; k<nz; k++){
          int idx = k * nw * nx * ny + j * nw * nx + i * nw + h;
          f(idx) = f4(w(h), x(i), y(j), z(k));
        }
      }
    }
  }
  Eigen::VectorXi ks(4), ns(4);
  ks << kw, kx, ky, kz;
  ns << nw, nx, ny, nz;
  Eigen::VectorXd wxyz(w.size() + x.size() + y.size() + z.size());
  wxyz << w, x, y, z;
  NdSpline sp = NdSpline(ks, ns, wxyz, f);

  nw = 10;
  nx = 10;
  ny = 10;
  nz = 10;
  w.resize(nw);
  x.resize(nx);
  y.resize(ny);
  z.resize(nz);
  f.resize(nw*nx*ny*nz);
  wxyz.resize(nw+nx+ny+nz);
  for (int i=0; i<nw; i++){
    w(i) = wmin + (double)(i) / (double)(nw-1) * (wmax-wmin);
  }
  for (int i=0; i<nx; i++){
    x(i) = xmin + (double)(i) / (double)(nx-1) * (xmax-xmin);
  }
  for (int i=0; i<ny; i++){
    y(i) = ymin + (double)(i) / (double)(ny-1) * (ymax-ymin);
  }
  for (int i=0; i<nz; i++){
    z(i) = zmin + (double)(i) / (double)(nz-1) * (zmax-zmin);
  }
  ns << nw, nx, ny, nz;
  wxyz << w, x, y, z;
  sp.Interpolate(ns, wxyz);

  double max_err = 0.0;
  std::ofstream outfile("file_4d.dat");
  outfile << "x, y, z, interpolated, true value" << std::endl;
  for(int h=0; h<nw; h++){
    for(int i=0; i<nx; i++){
      for(int j=0; j<ny; j++){
        for(int k=0; k<nz; k++){
          int idx = k * nw * nx * ny + j * nw * nx + i * nw + h;
          outfile
            << std::fixed << std::setw(8) << std::setprecision(4) << w(h)
            << std::fixed << std::setw(8) << std::setprecision(4) << x(i)
            << std::fixed << std::setw(8) << std::setprecision(4) << y(j)
            << std::fixed << std::setw(8) << std::setprecision(4) << z(k)
            << std::scientific << std::setw(18) << std::setprecision(6) << sp.fs(idx)
            << std::scientific << std::setw(18) << std::setprecision(6) << f4(w(h), x(i), y(j), z(k))
            << std::endl;
//          max_err = std::max(max_err, (std::abs(sp.fs(idx) - f4(w(h), x(i), y(j), z(k)))) / f4(w(h), x(i), y(j), z(k)));
          max_err = std::max(max_err, (std::abs(sp.fs(idx) - f4(w(h), x(i), y(j), z(k)))) );
        }
      }
    }
  }
  std::cout << "Maximum interpolation error: " << std::scientific << std::setw(18) << std::setprecision(6) << max_err << std::endl;
}

int main(int argc, char** argv)
{
  test_1d();
  test_2d();
  test_3d();
  test_4d();
}
