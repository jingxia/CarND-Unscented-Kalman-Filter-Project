#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  // calculate RMSE
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0; 

  // validate the input params
  if (estimations.size() != ground_truth.size() || estimations.size() == 0)
  {
    cout << "The input params are not valid." << endl;
    return rmse;
  }
  
  const size_t size = estimations.size();

  for (size_t i = 0; i != size; ++i)
  {
    VectorXd residuals = estimations[i] - ground_truth[i];
    residuals = residuals.array() * residuals.array();
    rmse += residuals;
  }

  rmse /= size;
  rmse = rmse.array().sqrt();

  return rmse;

}
