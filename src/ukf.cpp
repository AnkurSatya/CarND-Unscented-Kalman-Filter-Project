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
  std_a_ = 1.5;
  std_yawdd_ = M_PI/10.0;
  n_x_ = 5;
  n_aug_ = 7;

  P_ = MatrixXd(n_x_, n_x_);
  P_.fill(0.0);
  P_(0,0) = 0.5;//for px
  P_(1,1) = 0.5;//for py
  P_(2,2) = 1.0;//for v 
  P_(3,3) = 0.0;//for psi
  P_(4,4) = 0.0;//for psi'
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
      x_(2) = 0.0;//v
    }
    else
    {
      float p1 = meas_package.raw_measurements_(0);//rho
      float p2 = meas_package.raw_measurements_(1);//phi
      float p3 = meas_package.raw_measurements_(2);//rho'

      x_.fill(0.0);
      x_(0) = p1 * cos(p2),//px
      x_(1) = p1 * sin(p2),//py
      x_(2) = p3;//v
      }

    is_initialized_ = true;
    previous_timestamp = meas_package.timestamp_;
    return;
  }



  double delta_t = (meas_package.timestamp_ - previous_timestamp)/1000000.0;
  previous_timestamp = meas_package.timestamp_;
  
  //STEP 1 - Augmentation of the state vector and the state covariance matrix.
  VectorXd x_aug(n_aug_);
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_; 

  MatrixXd P_aug(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  //STEP 2 - Generating the sigma points.
  lambda_ = 3 - n_aug_;
  sigma_points = MatrixXd(n_aug_, 2*n_aug_+1);
  sigma_points.fill(0.0);

  sigma_points.col(0) = x_aug;
  MatrixXd sqrt_P_aug(n_aug_, n_aug_);
  sqrt_P_aug = P_aug.llt().matrixL();

  for(int i=0; i<n_aug_; i++)
  {
    // cout<<"Loop 1 - in"<<endl;
    sigma_points.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * sqrt_P_aug.col(i);
    sigma_points.col(i + n_aug_ + 1) = x_aug - sqrt(lambda_ + n_aug_) * sqrt_P_aug.col(i);
  }
 
  //STEP 3 - Prediction
  //Defining the sigma points weights.
  float value = 1/(2 * (lambda_ + n_aug_ )); 
  weights_ = VectorXd(2*n_aug_ + 1);

  weights_.fill(value);
  weights_(0) = lambda_/(lambda_ + n_aug_);

  Prediction(delta_t);
  use_laser_ = true;
  use_radar_ = true;

  // //STEP 4 - Update using the measurement values.
  if(use_laser_ == true && meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    UpdateLidar(meas_package);
  }
  else if(use_radar_== true && meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package);
  }
  cout<<"X -  "<<x_<<endl;
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
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);//necessary because otherwise it gets intialized with garbage values.
  
  for(int i=0; i<2*n_aug_ + 1; i++)
  {
    VectorXd col_val = sigma_points.col(i);
    VectorXd noises(n_x_);

    noises << 0.5 * delta_t*delta_t * cos(col_val(3)) * col_val(5),
              0.5 * delta_t*delta_t * sin(col_val(3)) * col_val(5),
              delta_t * col_val(5),
              0.5 * delta_t*delta_t * col_val(6),
              delta_t * col_val(6);
    Xsig_pred_.col(i) = col_val.head(n_x_) + noises;

    VectorXd change(n_x_);
    change.fill(0.0);

    if(fabs(col_val(4)) > 0.001)
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

  //STEP 2 - Predicting the mean and the covariance.
  x_.fill(0.0);

  for(int i=0; i<2*n_aug_ + 1; i++)
  {
    x_+= weights_(i) * Xsig_pred_.col(i);
  }

  P_.fill(0.0);
  for(int i=0; i<2*n_aug_ + 1; i++)
  {
    VectorXd diff = Xsig_pred_.col(i) - x_;
    tools.restrict_angle(diff(3));

    P_+= weights_(i) * diff * (diff.transpose());
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  int size = meas_package.raw_measurements_.size();
  MatrixXd z_sigma_pred(size, 2*n_aug_ + 1);
  z_sigma_pred.fill(0.0);

  //STEP 1 - Predicting the sigma points in the measurement space.
  for(int i=0; i<2*n_aug_ + 1; i++)
  {
    float px = sigma_points.col(i)(0);
    float py = sigma_points.col(i)(1);
    float vx = sigma_points.col(i)(2) * cos(sigma_points.col(i)(3));
    float vy = sigma_points.col(i)(2) * sin(sigma_points.col(i)(3));

    z_sigma_pred.col(i) << sqrt(px*px + py*py),
                           atan2(py, px),
                           (px*vx + py*vy)/sqrt(px*px + py*py);
  }

  //STEP 2 - Predicting the mean and the covariance.
  VectorXd z_(size);//Predicted state vector in the measurement space.
  z_.fill(0.0);

  for(int i=0; i<2*n_aug_ + 1; i++)
  {
    z_+=weights_(i) * z_sigma_pred.col(i);
  }

  MatrixXd S_(size, size);
  S_.fill(0.0);

  for(int i=0; i<2*n_aug_ + 1; i++)
  {
    VectorXd diff = z_sigma_pred.col(i) - z_;
    tools.restrict_angle(diff(1));

    S_+= weights_(i) * diff * (diff.transpose());
  }

  //STEP 3 - Defining the Measurement Noise Matrix and adding to the measurement noise variance - covariance matrix.
  MatrixXd R_(size, size);
  R_.fill(0.0);
  R_(0,0) = std_radr_ * std_radr_;
  R_(1,1) = std_radphi_ * std_radphi_;
  R_(2,2) = std_radrd_ * std_radrd_;

  S_+= R_;

  //STEP 4 - Evaluating the cross - correlation between sigma points in the measurement space and the state space and KALMAN GAIN.
  MatrixXd T_(n_x_, size);
  T_.fill(0.0);

  for(int i=0; i<2*n_aug_+1; i++)
  {
    VectorXd diff1 = Xsig_pred_.col(i) - x_;
    VectorXd diff2 = z_sigma_pred.col(i) - z_;
    tools.restrict_angle(diff2(1));

    T_+= weights_(i) * diff1 * (diff2.transpose()); 
  }

  MatrixXd K_(n_x_, size);
  K_.fill(0.0);

  K_ = T_ * (S_.inverse());//KALMAN GAIN

  //STEP 5 - Updating the state vector and the covariance - variance matrix.
  VectorXd diff_z = meas_package.raw_measurements_ - z_;
  tools.restrict_angle(diff_z(1)); 
  x_+= K_ * diff_z;
  P_-= K_ * S_ * (K_.transpose());

  //STEP 6 - Evaluating the NIS.
  double NIS;
  NIS = diff_z.transpose() * S_.inverse() * diff_z;
  NIS_list.push_back(NIS);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  int size = meas_package.raw_measurements_.size();

  MatrixXd H_(size, n_x_);
  H_.fill(0.0);
  H_(0,0) = 1;
  H_(1,1) = 1;

  VectorXd z_ = meas_package.raw_measurements_;
  VectorXd y_(size);
  MatrixXd R_(size, size);
  MatrixXd S_(size, size);
  MatrixXd K_(n_x_, size);
  MatrixXd I = MatrixXd::Identity(n_x_,n_x_);

  y_ = z_ - H_ * x_;
  
  R_.fill(0.0);
  R_(0,0) = std_laspx_ * std_laspx_;
  R_(1,1) = std_laspy_ * std_laspy_;

  S_ = H_ * P_ * (H_.transpose()) + R_;

  K_ = P_ * (H_.transpose()) * (S_.inverse());

  x_ = x_ + (K_ * y_);
  P_ = (I - (K_ * H_)) * P_;
}

vector<double> UKF::NIS_values()
{
  return NIS_list;
}