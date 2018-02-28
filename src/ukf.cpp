#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

namespace {
  const double EPSILON = 0.00001;
}

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.7;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  // state dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // state vector;
  x_ = VectorXd(n_x_);

  // state covariance matrix
  P_ = MatrixXd::Identity(n_x_, n_x_);
  P_(0,0) = 0.01;
  P_(1,1) = 0.01;
  

  // predicted sigma points
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // sigma point speading param
  lambda_ = 3 - n_x_;

  // weights vector
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_.fill(0.5 / (lambda_ + n_aug_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  // Augmented sigma points
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(const MeasurementPackage& meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  if (!is_initialized_)
  {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      double rho = meas_package.raw_measurements_[0];
      double fi = meas_package.raw_measurements_[1];
      double frdot = meas_package.raw_measurements_[2];
      
      x_ << rho * cos(fi), rho * sin(fi), frdot, 0, 0;
    }
    else 
    {
      double px = meas_package.raw_measurements_[0];
      double py = meas_package.raw_measurements_[1];

      x_ << px, py, 0, atan2(py, px), 0;
    }

    time_us_ = meas_package.timestamp_;

    is_initialized_ = true;

    return;
  }
  
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  // Predict
  if (fabs(dt) > EPSILON)
  {
    // Generate Sigma points
    GenerateAugmentedSigmaPoints();
    Prediction(dt);
  }

  // Update
  if (use_radar_  && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    UpdateRadar(meas_package);
  } 

  if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
    // Laser updates
    UpdateLidar(meas_package);
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  for (size_t i = 0; i != 2 * n_aug_ + 1; ++i)
  {
    VectorXd vec = VectorXd(n_x_);
    double px = Xsig_aug_(0, i);
    double py = Xsig_aug_(1, i);
    double v = Xsig_aug_(2, i);
    double fi = Xsig_aug_(3, i);
    double fi_dot = Xsig_aug_(4, i);
    double a = Xsig_aug_(5, i);
    double fi_dd = Xsig_aug_(6, i);

    if (fabs(fi_dot) < EPSILON)
    {
      vec(0) = px + v * cos(fi) * delta_t + a * pow(delta_t, 2) * cos(fi) / 2;
      vec(1) = py + v * sin(fi) * delta_t + a * pow(delta_t, 2) * sin(fi) / 2;
    }
    else
    {
      vec(0) = px + v * (sin(fi + fi_dot * delta_t) - sin(fi)) / fi_dot + a * pow(delta_t, 2) * cos(fi) / 2;
      vec(1) = py + v * (cos(fi) - cos(fi + fi_dot * delta_t)) / fi_dot + a * pow(delta_t, 2) * sin(fi) / 2;
    }

    vec(2) = v + a * delta_t;
    vec(3) = fi + fi_dot * delta_t + fi_dd * pow(delta_t, 2) / 2;
    vec(4) = fi_dot + fi_dd * delta_t;

    Xsig_pred_.col(i) = vec;
  }

  CalculatePredicatedMeanCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const MeasurementPackage& meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  const int n_z_ = 2;
  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  z_pred.fill(0.0);

  // measurements covariance matrix
  MatrixXd S = MatrixXd(n_z_, n_z_);
  S.fill(0.0);

  // predicted measurements matrix
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);

  for (size_t i = 0; i != 2 * n_aug_ + 1; ++i)
  {
    Zsig(0, i) = Xsig_pred_(0, i);
    Zsig(1, i) = Xsig_pred_(1, i);
  }

  // calculate mean
  for (size_t i = 0; i != 2 * n_aug_ + 1; ++i)
  {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  for (size_t i = 0; i != 2 * n_aug_ + 1; ++i)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * (z_diff * z_diff.transpose());
  }

  S(0, 0) += pow(std_laspx_, 2);
  S(1, 1) += pow(std_laspy_, 2);

  Update(meas_package.raw_measurements_, z_pred, S, Zsig);

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(const MeasurementPackage& meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  const int n_z_ = 3;
  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  z_pred.fill(0.0);

  // measurements covariance matrix
  MatrixXd S = MatrixXd(n_z_, n_z_);
  S.fill(0.0);

  // predicted measurements matrix
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);

  for (size_t i = 0; i != 2 * n_aug_ + 1; ++i)
  {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double fi = Xsig_pred_(3, i);
    double fi_dot = Xsig_pred_(4, i);

    double rou = sqrt(pow(px, 2) + pow(py, 2));
    double phi = atan2(py, px);
    double rou_dot = (px * cos(fi) * v + py * sin(fi) * v) / rou;

    Zsig(0, i) = rou;
    Zsig(1, i) = phi;
    Zsig(2, i) = rou_dot;
  }

  // calculate mean
  for (size_t i = 0; i != 2 * n_aug_ + 1; ++i)
  {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  NormalizeAngle(z_pred(1));

  for (size_t i = 0; i != 2 * n_aug_ + 1; ++i)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    NormalizeAngle(z_diff(1));

    S = S + weights_(i) * (z_diff * z_diff.transpose());
  }

  S(0, 0) += pow(std_radr_, 2);
  S(1, 1) += pow(std_radphi_, 2);
  S(2, 2) += pow(std_radrd_, 2);
  
  Update(meas_package.raw_measurements_, z_pred, S, Zsig);

}

void UKF::Update(const VectorXd& z, const VectorXd& z_pred, const MatrixXd& S, const MatrixXd& Zsig)
{

  const int n_z_ = z.size();

  MatrixXd Tc = MatrixXd(n_x_, n_z_);  
  Tc.fill(0.0);

  for (size_t i = 0; i != 2 * n_aug_ + 1; ++i)
  {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    NormalizeAngle(x_diff(3));

    VectorXd z_diff = Zsig.col(i) - z_pred;
    NormalizeAngle(z_diff(1));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  MatrixXd K = Tc * S.inverse();

  VectorXd z_diff = z - z_pred;
  NormalizeAngle(z_diff(1));

  x_ = x_ + K * z_diff;
  NormalizeAngle(x_(3));
  P_ = P_ - K * S * K.transpose();

}

/**
 * Calculate Augmented Sigma points
 */
void UKF::GenerateAugmentedSigmaPoints()
{
  // Augmented state vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0);

  // Augmented covariance matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0);

  Xsig_aug_.fill(0);

  x_aug.head(n_x_) = x_;

  P_aug.topLeftCorner(n_x_, n_x_) = P_;

  MatrixXd Q = MatrixXd(2, 2);
  Q << std_a_ * std_a_, 0,
       0, std_yawdd_ * std_yawdd_;
  
  P_aug.bottomRightCorner(2, 2) = Q;

  // create square root matrix
  MatrixXd sqt = P_aug.llt().matrixL();

  // fill in augmented sigma points
  Xsig_aug_.col(0) = x_aug;

  for (size_t i = 0; i != n_aug_; ++i)
  {
    Xsig_aug_.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * sqt.col(i);
    Xsig_aug_.col(i + n_aug_ + 1) = x_aug - sqrt(lambda_ + n_aug_) * sqt.col(i);
  }
}

void UKF::CalculatePredicatedMeanCovariance()
{
  x_.fill(0.0);
  for (size_t i = 0; i != 2 * n_aug_ + 1; ++i)
  {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  P_.fill(0.0);
  for (size_t i = 0; i != 2 * n_aug_ + 1; ++i)
  {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    NormalizeAngle(x_diff(3));

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}

void UKF::NormalizeAngle(double& angle)
{
    while (angle > M_PI)
    {
      angle -= 2 * M_PI;
    }
    while (angle < -1 * M_PI)
    {
      angle += 2 * M_PI;
    }
}
  
