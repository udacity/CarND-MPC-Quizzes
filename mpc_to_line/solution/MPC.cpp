#include "MPC.h"
#include <math.h>
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

using CppAD::AD;

// We set the number of timesteps to 25
// and the timestep evaluation frequency or evaluation
// period to 0.05.
size_t N = 25;
double dt = 0.05;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

double ref_v = 40;

// The solver takes all the state variables and actuator
// variables in a singular vector. Thus, we should to establish
// when one variable starts and another ends to make our lifes easier.
size_t delta_start = 0;
size_t a_start = delta_start + N;

typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

class FG_eval {
 public:
  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  Eigen::VectorXd coeffs;
  // Coefficients of the fitted polynomial.
  FG_eval(Eigen::VectorXd coeffs, double x0, double y0, double psi0, double v0) {
    this->coeffs = coeffs;
    this->x0 = x0;
    this->y0 = y0;
    this->psi0 = psi0;
    this->v0 = v0;
  }
  
  double x0, y0, psi0, v0;

  void run_car_model(const ADvector& vars, ADvector& x_out, ADvector& y_out, ADvector& psi_out, ADvector& v_out)
  {
    AD<double> x = x0, y = y0, psi = psi0, v = v0;
    for(int t=0; t< N; t++){
      auto delta = vars[delta_start + t];
      auto a = vars[a_start + t];

      x = x + v * CppAD::cos(psi) * dt;
      y = y + v * CppAD::sin(psi) * dt;
      psi = psi + v * delta / Lf * dt;
      v = v + a * dt;

      x_out[t] = x;
      y_out[t] = y;
      psi_out[t] = psi;
      v_out[t] = v;

    }
  }

  // `fg` is a vector containing the cost and constraints.
  // `vars` is a vector containing the variable values (actuators).
  void operator()(ADvector& fg, const ADvector& vars) {
    ADvector x_pred(N), y_pred(N), psi_pred(N), v_pred(N);
    run_car_model(vars, x_pred, y_pred, psi_pred, v_pred);

    // The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
    fg[0] = 0;

    // The part of the cost based on the reference state.
    for (int t = 0; t < N; t++) {
      // CTE
      AD<double> f0 = coeffs[0] + coeffs[1] * x_pred[t];
      fg[0] += CppAD::pow(f0 - y_pred[t], 2);
      // epsi
      AD<double> psides0 = CppAD::atan(coeffs[1]);
      fg[0] += CppAD::pow(psides0 - psi_pred[t], 2);

      fg[0] += CppAD::pow(ref_v - v_pred[t], 2);
    }

    // Minimize the use of actuators.
    for (int t = 0; t < N ; t++) {
      fg[0] += CppAD::pow(vars[delta_start + t], 2);
      fg[0] += CppAD::pow(vars[a_start + t], 2);
    }

    // Minimize the value gap between sequential actuations.
    for (int t = 0; t < N - 1; t++) {
      fg[0] += CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
      fg[0] += CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
    }

    //
    // Setup Constraints
    //
    
    // NONE !!!
  }
};

//
// MPC class definition
//

MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd x0, Eigen::VectorXd coeffs) {
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  double x = x0[0];
  double y = x0[1];
  double psi = x0[2];
  double v = x0[3];

  // number of independent variables
  // N timesteps == N - 1 actuations
  size_t n_vars = (N) * 2;
  // Number of constraints
  size_t n_constraints = 0;

  // Initial value of the independent variables.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0.0;
  }

  // Lower and upper limits for x
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  // NOTE: Feel free to change this to something else.
  for (int i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }

  // Acceleration/decceleration upper and lower limits.
  // NOTE: Feel free to change this to something else.
  for (int i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for constraints
  // All of these should be 0 except the initial
  // state indices.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) { // XXX loop of length 0 
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  // Object that computes objective and constraints
  FG_eval fg_eval(coeffs, x, y, psi, v);

  // options
  std::string options;
  options += "Integer print_level  0\n";
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  //
  // Check some of the solution values
  //
  bool ok = true;
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  FG_eval::ADvector vars_ster(n_vars), x_pred(N), y_pred(N), psi_pred(N), v_pred(N);
  for(int i=0; i<n_vars; i++) vars_ster[i] = solution.x[i];
  fg_eval.run_car_model(vars_ster, x_pred, y_pred, psi_pred, v_pred);
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;
  return {CppAD::Value(x_pred[0]), CppAD::Value(y_pred[0]),
    CppAD::Value(psi_pred[0]), CppAD::Value(v_pred[0]), 
    solution.x[delta_start], solution.x[a_start]};
}

//
// Helper functions to fit and evaluate polynomials.
//

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

int main() {
  MPC mpc;
  int iters = 50;

  Eigen::VectorXd ptsx(2);
  Eigen::VectorXd ptsy(2);
  ptsx << -100, 100;
  ptsy << -1, -1;

  // The polynomial is fitted to a straight line so a polynomial with
  // order 1 is sufficient.
  auto coeffs = polyfit(ptsx, ptsy, 1);

  // NOTE: free feel to play around with these
  double x = -1;
  double y = 10;
  double psi = 0;
  double v = 10;
  // The cross track error is calculated by evaluating at polynomial at x, f(x)
  // and subtracting y.
  double cte = polyeval(coeffs, x) - y;
  // Due to the sign starting at 0, the orientation error is -f'(x).
  // derivative of coeffs[0] + coeffs[1] * x -> coeffs[1]
  double epsi = psi - atan(coeffs[1]);

  Eigen::VectorXd state(4);
  state << x, y, psi, v;

  std::vector<double> x_vals = {state[0]};
  std::vector<double> y_vals = {state[1]};
  std::vector<double> psi_vals = {state[2]};
  std::vector<double> v_vals = {state[3]};
  std::vector<double> cte_vals = {cte};
  std::vector<double> epsi_vals = {epsi};
  std::vector<double> delta_vals = {};
  std::vector<double> a_vals = {};

  for (size_t i = 0; i < iters; i++) {
    std::cout << "Iteration " << i << std::endl;

    auto vars = mpc.Solve(state, coeffs);

    x_vals.push_back(vars[0]);
    y_vals.push_back(vars[1]);
    psi_vals.push_back(vars[2]);
    v_vals.push_back(vars[3]);
  
    double cte = polyeval(coeffs, vars[0]) - vars[1];
    double epsi = vars[2] - atan(coeffs[1]);


    cte_vals.push_back(cte);
    epsi_vals.push_back(epsi);

    delta_vals.push_back(vars[4]);
    a_vals.push_back(vars[5]);

    state << vars[0], vars[1], vars[2], vars[3];
    std::cout << "x = " << vars[0] << std::endl;
    std::cout << "y = " << vars[1] << std::endl;
    std::cout << "psi = " << vars[2] << std::endl;
    std::cout << "v = " << vars[3] << std::endl;
    std::cout << "cte = " << cte << std::endl;
    std::cout << "epsi = " << epsi << std::endl;
    std::cout << "delta = " << vars[4] << std::endl;
    std::cout << "a = " << vars[5] << std::endl;
    std::cout << std::endl;
  }

  // Plot values
  // NOTE: feel free to play around with this.
  // It's useful for debugging!
  plt::subplot(3, 1, 1);
  plt::title("CTE");
  plt::plot(cte_vals);
  plt::subplot(3, 1, 2);
  plt::title("Delta (Radians)");
  plt::plot(delta_vals);
  plt::subplot(3, 1, 3);
  plt::title("Velocity");
  plt::plot(v_vals);

  plt::show();
}
