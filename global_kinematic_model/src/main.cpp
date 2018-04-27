// In this quiz you'll implement the global kinematic model.
#include <math.h>
#include <iostream>
#include "Eigen-3.3/Eigen/Core"
#include <sstream>
#include <fstream>

//
// Helper functions
//
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

const double Lf = 2;

// TODO: Implement the global kinematic model.
// Return the next state.
//
// NOTE: state is [x, y, psi, v]
// NOTE: actuators is [delta, a]
Eigen::VectorXd globalKinematic(Eigen::VectorXd state,
                                Eigen::VectorXd actuators, double dt, double Lf) {

  return next_state;
}

void write_log(std::string namefile, int n, Eigen::VectorXd state,Eigen::VectorXd actuators, double dt, double Lf) {

  std::ofstream dataFile;
  dataFile.open(namefile);

  dataFile << state[0] << " " << state[1] << " " <<  state[2]<< " " <<  state[3];
  dataFile << "\n";
  for(int i=0; i<n;i++){
      Eigen::VectorXd next_state = globalKinematic(state, actuators, dt, Lf);
      state << next_state[0], next_state[1], next_state[2],next_state[3];
      dataFile << next_state[0] << " " << next_state[1] << " " <<  next_state[2]<< " " <<  next_state[3];
      dataFile << "\n";
      
  }
  
}

int main() {
  // [x, y, psi, v]
  Eigen::VectorXd state(4);
  // [delta, v]
  Eigen::VectorXd actuators(2);

  state << 0, 0, deg2rad(0), 0;
  actuators << deg2rad(5), 1;

  std::string path = "../logs/"; 
  std::string namefile = path + "Lf1.txt";
  write_log(namefile, 75, state, actuators, 0.3, 1);
  namefile = path + "Lf2.txt";
  write_log(namefile, 75, state, actuators, 0.3, 2.67);
  namefile = path + "Lf4.txt";
  write_log(namefile, 75, state, actuators, 0.3, 4);
}