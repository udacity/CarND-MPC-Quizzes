// In this quiz you'll implement the dynamic constraints
// for the state [cte, epsi, v].
//
// The cost is already setup and the initial state is already
// constrained. You may play around with the cost :-)
#include "MPC.h"

#include <cmath>
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include <vector>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

namespace {
using CppAD::AD;

size_t N = 25;
double dt = 0.10;
double Lf = 2.67;

double ref_cte = 0;
double ref_epsi = 0;
// Feel free to change this.
double ref_v = 30;

size_t cte_start = 0;
size_t epsi_start = cte_start + N;
size_t v_start = epsi_start + N;
size_t sa_start = v_start + N;
size_t a_start = sa_start + N - 1;

class FG_eval {
 public:
  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& x) {
    // Total Cost
    fg[0] = 0;

    // Control Cost
    for (int i = 0; i < N - 1; i++) {
      fg[0] += CppAD::pow(x[sa_start + i], 2);
      fg[0] += CppAD::pow(x[a_start + i], 2);
    }

    // Deriv Input Cost
    for (int i = 0; i < N - 2; i++) {
      fg[0] += CppAD::pow(x[sa_start + i + 1] - x[sa_start + i], 2);
      fg[0] += CppAD::pow(x[a_start + i + 1] - x[a_start + i], 2);
    }

    // State Cost
    for (int i = 0; i < N; i++) {
      fg[0] += CppAD::pow(x[cte_start + i] - ref_cte, 2);
      fg[0] += CppAD::pow(x[epsi_start + i] - ref_epsi, 2);
      fg[0] += CppAD::pow(x[v_start + i] - ref_v, 2);
    }

    // 0 - cost
    // Initial state constraints.
    fg[1] = x[cte_start];
    fg[N + 1] = x[epsi_start];
    fg[2 * N + 1] = x[v_start];

    // NOTE 1: Use CppAD::sin instead of regular C++ sin.
    // NOTE 2: Be mindful of the above initial constraints indexes while
    // you're setting up the below dynamic constraint.
    //
    // TODO: Setup constraints
    //
    // cte[t+1] - (cte[t] + v[t] * sin(epsi[t] + sa[t]) * dt) == 0
    // epsi[t+1] - (epsi[t] + v[t] * sin(sa[t]) * dt) == 0
    // v[t+1] - (v[t] + a[t] * dt) == 0
  }
};
}

//
// MPC class definition
//

MPC::MPC() {}
MPC::~MPC() {}

tuple<vector<double>, vector<double>, double> MPC::Solve(vector<double> x0) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  double cte = x0[0];
  double epsi = x0[1];
  double v = x0[2];

  // number of independent variables (domain dimension for f and g)
  size_t nx = N * 3 + (N - 1) * 2;
  // number of constraints (range dimension for g)
  size_t ng = N * 3;

  // initial value of the independent variables
  Dvector xi(nx);
  for (int i = 0; i < nx; i++) {
    xi[i] = 0.0;
  }
  xi[0] = cte;
  xi[25] = epsi;
  xi[50] = v;

  // lower and upper limits for the state

  Dvector xl(nx), xu(nx);
  for (int i = 0; i < sa_start; i++) {
    xl[i] = -1.0e19;
    xu[i] = 1.0e19;
  }
  // [-25, 25] degrees steering angle
  // Change this if you like.
  for (int i = sa_start; i < a_start; i++) {
    xl[i] = -0.436332;
    xu[i] = 0.436332;
  }
  for (int i = a_start; i < nx; i++) {
    xl[i] = -1.0;
    xu[i] = 1.0;
  }

  // lower and upper limits for constraints
  Dvector gl(ng), gu(ng);
  for (int i = 0; i < ng; i++) {
    gl[i] = 0;
    gu[i] = 0;
  }

  // initial state constraint
  gl[0] = cte;
  gu[0] = cte;
  gl[25] = epsi;
  gu[25] = epsi;
  gl[50] = v;
  gu[50] = v;

  // object that computes objective and constraints
  FG_eval fg_eval;

  // options
  std::string options;
  // turn off any printing
  options += "Integer print_level 0\n";
  options += "Sparse  true  forward\n";
  options += "Sparse  true  reverse\n";
  options += "Numeric max_cpu_time  0.05\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(options, xi, xl, xu, gl, gu, fg_eval,
                                        solution);
  //
  // Check some of the solution values
  //
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  std::cout << "CTE " << solution.x[cte_start + 1] << std::endl;
  std::cout << "Epsi " << solution.x[epsi_start + 1] << std::endl;
  std::cout << "Velocity " << solution.x[v_start + 1] << std::endl;
  std::cout << "Steering Angle " << solution.x[sa_start] << std::endl;
  std::cout << "Acceleration " << solution.x[a_start] << std::endl;
  std::cout << std::endl;

  auto x1 = {solution.x[cte_start + 1], solution.x[epsi_start + 1],
             solution.x[v_start + 1]};
  auto u0 = {solution.x[sa_start], solution.x[a_start]};
  auto cost = solution.obj_value;
  return std::make_tuple(x1, u0, cost);
}

int main() {
  MPC mpc;
  auto iters = 50;
  std::vector<double> x = {5., 0., 0.};
  std::vector<double> u = {0., 0.};

  std::vector<double> cte_vals = {x[0]};
  std::vector<double> epsi_vals = {x[1]};
  std::vector<double> v_vals = {x[2]};
  std::vector<double> sa_vals = {u[0]};
  std::vector<double> a_vals = {u[1]};

  auto vals = mpc.Solve(x);

  for (size_t i = 0; i < iters; i++) {
    std::cout << "Iteration " << i << std::endl;

    vals = mpc.Solve(x);
    x = std::get<0>(vals);
    u = std::get<1>(vals);

    cte_vals.push_back(x[0]);
    epsi_vals.push_back(x[1]);
    v_vals.push_back(x[2]);
    sa_vals.push_back(u[0]);
    a_vals.push_back(u[1]);
  }

  // Plot values
  plt::subplot(3, 1, 1);
  plt::title("CTE");
  plt::plot(cte_vals);
  plt::subplot(3, 1, 2);
  plt::title("Epsi (Radians)");
  plt::plot(sa_vals);
  plt::subplot(3, 1, 3);
  plt::title("Velocity");
  plt::plot(v_vals);

  plt::show();
}