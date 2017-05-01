// In this quiz you'll implement the global kinematic model.

#include <math.h>
#include <iostream>
#include "Eigen-3.3/Eigen/Core"

constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

const double Lf = 2;

// TODO: Implement the global kinematic model.
// Return the next state.
//
// NOTE: state is [x, y, psi, v]
// NOTE: actuators is [sa, a]
Eigen::VectorXd globalKinematic(Eigen::VectorXd state,
                                Eigen::VectorXd actuators, double dt) {
  Eigen::VectorXd next_state(state.size());

  auto x = state(0);
  auto y = state(1);
  auto psi = state(2);
  auto v = state(3);

  auto sa = actuators(0);
  auto a = actuators(1);

  next_state(0) = x + v * cos(psi) * dt;
  next_state(1) = y + v * sin(psi) * dt;
  next_state(2) = psi + v / Lf * sa * dt;
  next_state(3) = v + a * dt;
  return next_state;
}

int main() {
  Eigen::VectorXd state(4);
  Eigen::VectorXd actuators(2);
  state << 0, 0, deg2rad(45), 1;
  actuators << deg2rad(5), 1;

  auto next_state = globalKinematic(state, actuators, 0.3);

  std::cout << next_state << std::endl;
}