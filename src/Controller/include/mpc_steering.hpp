#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/QR>
#include <eigen3/Eigen/Householder>
#include <iostream>
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
    #include "common.hpp"

#define PLOT (1)
using namespace std;
namespace plt = matplotlibcpp;

// namespace AD = CppAD;
class MPC{
public:
    MPC(){};
    ~MPC(){};
    std::vector<double> mpc_solve(Eigen::VectorXd state,Eigen::VectorXd coeff,
                                  std::vector<std::vector<double>>& obstacles);
    
};


