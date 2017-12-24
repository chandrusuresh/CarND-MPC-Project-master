#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC {
 public:
  MPC();

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
    string outFileName;
    int inputDelayPts;
    double prev_delta;
    double prev_a;
    double latency;
    void writeToFile(Eigen::VectorXd state, double px, double py, double steer, double thr);
};

#endif /* MPC_H */
