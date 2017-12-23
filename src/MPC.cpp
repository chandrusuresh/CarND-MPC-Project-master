#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

#include <math.h>
#include "Eigen-3.3/Eigen/QR"

using CppAD::AD;

// TODO: Set the timestep length and duration
size_t N = 50;
double dt = 0.01;

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
double v_ref = 30;
size_t x_start = 0;
size_t y_start = x_start + N;
size_t psi_start = y_start + N;
size_t v_start = psi_start + N;
size_t cte_start = v_start + N;
size_t epsi_start = cte_start + N;
size_t delta_start = epsi_start + N;
size_t a_start = delta_start + N-1;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
      
      // Cost Function
      fg[0] = 0.0;
      
      // State Cost
      for (int t=0; t < N; t++)
      {
          fg[0] += CppAD::pow(vars[cte_start+t],2);
          fg[0] += CppAD::pow(vars[epsi_start+t],2);
          fg[0] += CppAD::pow(vars[v_start+t]-v_ref,2);
      }
      
      // Input cost
      for (int t=0; t < N-1; t++)
      {
          fg[0] += 100*CppAD::pow(vars[delta_start+t],2);
          fg[0] += CppAD::pow(vars[a_start+t],2);
      }
      
      for (int t=1; t < N-1; t++)
      {
          fg[0] += 10*CppAD::pow(vars[delta_start+t] - vars[delta_start+t-1],2);
          fg[0] += CppAD::pow(vars[a_start+t] - vars[a_start+t-1],2);
      }
      
//      // Cost Function
//      fg[0] = 0.0;
//
//      // State Cost
//      for (int t=0; t < N; t++)
//      {
//          fg[0] += CppAD::pow(vars[cte_start+t],2);
//          fg[0] += 100*CppAD::pow(vars[epsi_start+t],2);
//          fg[0] += CppAD::pow(vars[v_start+t]-v_ref,2);
//      }
//
//      // Input cost
//      for (int t=0; t < N-1; t++)
//      {
//          fg[0] += 100*CppAD::pow(vars[delta_start+t],2);
//          fg[0] += CppAD::pow(vars[a_start+t],2);
//      }
//
//      for (int t=1; t < N-1; t++)
//      {
//          fg[0] += 10*CppAD::pow(vars[delta_start+t] - vars[delta_start+t-1],2);
//          fg[0] += 10*CppAD::pow(vars[a_start+t] - vars[a_start+t-1],2);
//      }
      
      // Constraints
      fg[1+x_start] = vars[x_start];
      fg[1+y_start] = vars[y_start];
      fg[1+psi_start] = vars[psi_start];
      fg[1+v_start] = vars[v_start];
      fg[1+cte_start] = vars[cte_start];
      fg[1+epsi_start] = vars[epsi_start];
      
      for (int t=1; t < N; t++)
      {
          AD<double> x0 = vars[x_start+t-1];
          AD<double> y0 = vars[y_start+t-1];
          AD<double> psi0 = vars[psi_start+t-1];
          AD<double> v0 = vars[v_start+t-1];
          AD<double> cte0 = vars[cte_start+t-1];
          AD<double> epsi0 = vars[epsi_start+t-1];
          AD<double> delta0 = vars[delta_start+t-1];
          AD<double> a0 = vars[a_start+t-1];
          
          AD<double> x1 = vars[x_start+t];
          AD<double> y1 = vars[y_start+t];
          AD<double> psi1 = vars[psi_start+t];
          AD<double> v1 = vars[v_start+t];
          AD<double> cte1 = vars[cte_start+t];
          AD<double> epsi1 = vars[epsi_start+t];
          
          AD<double> f0 = 0;
          AD<double> tan_psides0 = 0;
          for (int i=0; i < coeffs.size();i++)
          {
              f0 += coeffs[i]*CppAD::pow(x0,i);
              if (i > 0)
              {
                  tan_psides0 += i*CppAD::pow(x0,i-1)*coeffs[i];
              }
          }
          AD<double> psides0 = CppAD::atan(tan_psides0);
          
          fg[1+x_start+t]    =   x1  - (x0 + v0 * CppAD::cos(psi0) * dt);
          fg[1+y_start+t]    =   y1  - (y0 + v0 * CppAD::sin(psi0) * dt);
          fg[1+psi_start+t]  =  psi1 - (psi0 + v0 * delta0/Lf* dt);
          fg[1+v_start+t]    =   v1  - (v0 + a0 * dt);
          fg[1+cte_start+t]  =  cte1 - ((f0 - y0) + v0 * CppAD::sin(psi0) * dt);
          fg[1+epsi_start+t] = epsi1 - ((psides0 - epsi0) + v0 * delta0/Lf* dt);
      }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {
    this->inputDelayPts = 15;
    this->outFileName = "MPC_refVel_" + std::to_string(v_ref) + "_inputDelayed_" + std::to_string(this->inputDelayPts) + "_latency_thrTuning.csv";
    ofstream myfile;
    myfile.open(outFileName);//,ios::app);
    myfile << "x,y,psi,v,cte,epsi,steer,throttle" << std::endl;
}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

    double x = state[0];
    double y = state[1];
    double psi = state[2];
    double v = state[3];
    double cte = state[4];
    double epsi = state[5];
    
  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  size_t n_vars = 6*N + 2*(N-1);
  // TODO: Set the number of constraints
  size_t n_constraints = 6*N;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }
    vars[x_start] = x;
    vars[y_start] = y;
    vars[psi_start] = psi;
    vars[v_start] = v;
    vars[cte_start] = cte;
    vars[epsi_start] = epsi;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.
    for (i = 0; i < n_vars; i++) {
        vars_lowerbound[i] = -1.0e19;
        vars_upperbound[i] =  1.0e19;
    }
    
    for (i = delta_start; i < a_start; i++) {
        vars_lowerbound[i] = -0.43633;
        vars_upperbound[i] =  0.43633;
    }
    for (i = a_start; i < n_vars; i++) {
        vars_lowerbound[i] = -1.0;
        vars_upperbound[i] =  1.0;
    }
    
  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0.0;
    constraints_upperbound[i] = 0.0;
  }
    constraints_lowerbound[x_start] = x;
    constraints_lowerbound[y_start] = y;
    constraints_lowerbound[psi_start] = psi;
    constraints_lowerbound[v_start] = v;
    constraints_lowerbound[cte_start] = cte;
    constraints_lowerbound[epsi_start] = epsi;
    
    constraints_upperbound[x_start] = x;
    constraints_upperbound[y_start] = y;
    constraints_upperbound[psi_start] = psi;
    constraints_upperbound[v_start] = v;
    constraints_upperbound[cte_start] = cte;
    constraints_upperbound[epsi_start] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
    vector<double> sol;
    for (i=x_start; i < x_start+N; i++)
    {
        sol.push_back(solution.x[i]);
    }
    for (i=y_start; i < y_start+N; i++)
    {
        sol.push_back(solution.x[i]);
    }
    for (i=delta_start; i < delta_start+N-1; i++)
    {
        sol.push_back(solution.x[i]);
    }
    for (i=a_start; i < a_start+N-1; i++)
    {
        sol.push_back(solution.x[i]);
    }
    return {solution.x[x_start + 1], solution.x[x_start + 1 + 1*N/5], solution.x[x_start + 1 + 2*N/5], solution.x[x_start + 1 + 3*N/5], solution.x[x_start + 1 + 4*N/5],
        solution.x[y_start + 1], solution.x[y_start + 1 + 1*N/5], solution.x[y_start + 1 + 2*N/5], solution.x[y_start + 1 + 3*N/5], solution.x[y_start + 1 + 4*N/5],
        solution.x[delta_start + this->inputDelayPts],   solution.x[a_start + this->inputDelayPts]};
}

void MPC::writeToFile(Eigen::VectorXd state, double px, double py, double steer, double thr) {
    ofstream myfile;
    myfile.open(this->outFileName,ios::app);
    myfile << px << "," << py << ",";
    for (int i=2; i < state.size(); i++)
    {
        myfile << state[i] <<  ",";
    }
    myfile << steer << "," << thr << std::endl;
    myfile.close();
}

