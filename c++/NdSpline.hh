#ifndef NdSpline_hh
#define NdSpline_hh 1
#include <iostream>
#include <vector>
#include <Eigen/Dense>

class NdSpline
{
  public:
    ~NdSpline();
    NdSpline(Eigen::VectorXi ks, Eigen::VectorXi ns, Eigen::VectorXd x, Eigen::VectorXd f);
    void Interpolate(const Eigen::VectorXi&, const Eigen::VectorXd&);
    Eigen::VectorXd fs;

  private:
    class abscissa {
      public:
        Eigen::VectorXd x, t;
        int nx = 0;
        int nt = 0;
        int k = 0;
        bool extrapolation = true;

        ~abscissa(){};
        abscissa(int, const Eigen::VectorXd&);
        void set_knot_vector();
        int find_interval(double) const;
        double basis_function(int, int, double) const;
        double basis_function_rec(int, int, int, double) const;
    };

    class abscissa_interpolant {
      public:
        Eigen::VectorXd x;
        int n_interpolant = 0;

        ~abscissa_interpolant(){};
        abscissa_interpolant(const Eigen::VectorXd&);
    };

    class grid_interpolant {
      public:
        std::vector<abscissa_interpolant> Ndrctn;
        int ndim = 0;
        int n_interpolant = 0;

        ~grid_interpolant(){};
        grid_interpolant(const Eigen::VectorXi&, const Eigen::VectorXd&);
    };

    std::vector<abscissa> Ndrctn;
    Eigen::VectorXd coefs;
    int ndim = 0;
    int nf_in = 0;

    void set_coefs();
    void interpolate_1d(const abscissa&, const Eigen::VectorXd&, const Eigen::MatrixXd&, Eigen::MatrixXd&);
    void set_coefs_1d(const abscissa&, const Eigen::MatrixXd&, Eigen::MatrixXd&);

    //void erase_abscissa(abscissa&);
    //void initialize_abscissa(abscissa&, int, Eigen::VectorXd&);
    //void set_knot_vector(abscissa&);
    //int find_interval(const abscissa&, const Eigen::VectorXd&);
    //double basis_function(const abscissa&, int, int, const Eigen::VectorXd&);

    void erase_grid_interpolant(grid_interpolant&);
    void initialize_grid_interpolant(grid_interpolant&, const Eigen::VectorXd&, const Eigen::VectorXd&);

    void erase_abscissa_interpolant(abscissa_interpolant&);
    void initialize_abscissa_interpolant(abscissa_interpolant&, const Eigen::VectorXd&);
};
#endif
