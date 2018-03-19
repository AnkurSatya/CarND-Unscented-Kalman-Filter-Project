#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
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
  is_initialized = false;

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;
  previous_timestamp = 0;
  std_a_ = 6.0;
  std_yawdd_ = M_PI/8.0;
  n_x_ = 5;
  n_aug_ = 7;

  P_ = MatrixXd(n_x_, n_x_);
  P_.fill(0.0);
  P_(0,0) = 1.0;//for px
  P_(1,1) = 1.0;//for py
  P_(2,2) = 25.0;//for v //std_dev for v = 5
  P_(3,3) = 1.0;//for psi
  P_(4,4) = 1.0;//for psi'
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  //STEP 0 - Initializing the vectors and matrices.
  if(!is_initialized_)
  {
    if(meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      x_.fill(0.0);
      x_.head(2) = meas_package.raw_measurements_;//px, py
      x_(2) = 1.0;//v
    }
    else
    {
      float p1 = meas_package.raw_measurements_(0);//rho
      float p2 = meas_package.raw_measurements_(1);//phi
      float p3 = meas_package.raw_measurements_(2);//rho'

      x_.fill(0.0);
      x_ << p1 * cos(p2),//px
            p1 * sin(p2),//py
            p3;//v
      }

    is_initialized_ = true;
    previous_timestamp = meas_package.timestamp_;
    cout<<"X - "<<x_<<endl;
    cout<<"P - "<<P_<<endl;
    return;
  }



  double delta_t = (meas_package.timestamp_ - previous_timestamp)/1000000.0;
  previous_timestamp = meas_package.timestamp_;
  
  //STEP 1 - Augmentation of the state vector and the state covariance matrix.
  VectorXd x_aug(n_aug_);
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_.head(n_x_); 

  MatrixXd P_aug(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  //STEP 2 - Generating the sigma points.
  lambda_ = 3 - n_aug_;
  sigma_points = MatrixXd(n_aug_, 2*n_aug_+1);

  sigma_points.col(0) = x_aug;
  MatrixXd sqrt_P_aug = P_aug.llt().matrixL();

  for(int i=0; i<n_aug_; i++)
  {
    cout<<"Loop 1 - in"<<endl;
    sigma_points.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * sqrt_P_aug.col(i);
    sigma_points.col(i + n_aug_ + 1) = x_aug - sqrt(lambda_ + n_aug_) * sqrt_P_aug.col(i);
  }
  cout<<"Loop 1 - out"<<endl;

  //STEP 3 - Prediction
  //Defining the sigma points weights.
  float value = 1/(2 * (lambda_ + 2 * n_aug_ + 1)); 
  weights_ = VectorXd(2*n_aug_ + 1);

  weights_.fill(value);
  weights_(0) = lambda_/(lambda_ + 2 * n_aug_ + 1);

  // Prediction(delta_t);

  //STEP 4 - Update using the measurement values,
  cout<<"X - "<<x_<<endl;
  cout<<"P - "<<P_<<endl;
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

  //STEP 1 - Predicting the sigma points.
  for(int i=0; i<2*n_aug_ + 1; i++)
  {
    cout<<"Loop 2 - in"<<endl;
    VectorXd col_val = sigma_points.col(i);
    VectorXd noises(n_x_);
    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

    noises << 0.5 * delta_t*delta_t * cos(col_val(3)) * col_val(5),
              0.5 * delta_t*delta_t * sin(col_val(3)) * col_val(5),
              delta_t * col_val(5),
              0.5 * delta_t*delta_t * col_val(6),
              delta_t * col_val(6);
    Xsig_pred_.col(i) = col_val.head(n_x_) + noises;

    VectorXd change(n_x_);

    if(col_val(4) != 0)
    {
      change << (col_val(2)/col_val(4)) * (sin(col_val(3) + col_val(4) * delta_t) - sin(col_val(3))),
                (col_val(2)/col_val(4)) * (-cos(col_val(3) + col_val(4) * delta_t) + cos(col_val(3))),
                0,
                col_val(4) * delta_t,
                0;
    }
    else
    {
      change << col_val(2) * cos(col_val(3)) * delta_t,
                col_val(2) * sin(col_val(3)) * delta_t,
                0,
                col_val(4) * delta_t,
                0;
    }
    Xsig_pred_.col(i)+= change;
  }
  cout<<"Loop 2 - out"<<endl;

  //STEP 2 - Predicting the mean and the covariance.
  x_.fill(0.0);

  for(int i=0; i<2*n_aug_ + 1; i++)
  {
    cout<<"Loop 3 - in"<<endl;
    x_+= weights_(i) * Xsig_pred_.col(i);
  }
  cout<<"Loop 3 - out"<<endl;

  P_.fill(0.0);
  for(int i=0; i<2*n_aug_ + 1; i++)
  {
    cout<<"Loop 4 - in"<<endl;
    P_+= weights_(i) * (Xsig_pred_.col(i) - x_) * ((Xsig_pred_.col(i) - x_).transpose());
  }
  cout<<"Loop 4 - out"<<endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
