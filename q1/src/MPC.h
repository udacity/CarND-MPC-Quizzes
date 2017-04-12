#ifndef MPC_H
#define MPC_H

#include <vector>

using namespace std;

class MPC {
 public:
  MPC();

  virtual ~MPC();

  // Solve the model given an initial state.
  // Return the next state, inputs and cost.
  tuple<vector<double>, vector<double>, double> Solve(vector<double> x0);
};

#endif /* MPC_H */
