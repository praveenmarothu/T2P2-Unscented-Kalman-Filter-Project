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
  std_a_ = 1.60;
  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.60;
  
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
    
  //set state dimension
  n_x_ = x_.size();
    
  //set augmented dimension
  n_aug_ = n_x_ + 2;
    
  //define spreading parameter
  lambda_ = 3 - n_aug_;
  
  //set weights dimensions
  weights_ = VectorXd( 2*n_aug_ + 1 );
    
  //store previous timestamp
  time_us_ = 0;
  
  //sigma point matrix
  Xsig_pred_  = MatrixXd(n_x_, 2*n_aug_ + 1);
  
  is_initialized_ = false;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

    if(!is_initialized_){
       
        //Init covariance matrix
        P_ << 1.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 1.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 1.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 1.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 1.0;
        
        
        if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){

            float ro = meas_package.raw_measurements_[0];
            float phi = meas_package.raw_measurements_[1];            
            float ro_dot = meas_package.raw_measurements_[2];
            
            //polar to cartesian
            float p_x = ro * cos(phi);
            float p_y = ro * sin(phi);
            float vx = ro_dot * cos(phi);
            float vy = ro_dot * sin(phi);
            float v = sqrt(vx*vx + vy*vy);

            //Initial location and velocity
            x_ << p_x, p_y, v, 0, 0;
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
        {
            
            float p_x = meas_package.raw_measurements_[0];
            float p_y = meas_package.raw_measurements_[1];
            float v  = 0;

            //Initial location and zero velocity
            x_ << p_x, p_y, v, 0, 0;

        }       

       weights_(0) = lambda_/(lambda_+n_aug_);
       for (int i=1; i< weights_.size(); i++) 
       {
            weights_(i) = 0.5/(lambda_+n_aug_);
       }

        time_us_ = meas_package.timestamp_;
        is_initialized_ = true;
        return;
    }


    //----------------------------------------------------------------------------
    //   PREDICTION
    //----------------------------------------------------------------------------
    
    //delta t
    float delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
    time_us_ = meas_package.timestamp_;

    if( delta_t > 0.001) Prediction(delta_t);
    
    
    if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) UpdateLidar(meas_package);
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) UpdateRadar(meas_package);
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */


void UKF::Prediction(double delta_t){

    //----------------------------------------------------------------------------
    //   CREATE SIGMA POINTS
    //----------------------------------------------------------------------------

    //augmented mean vector
    VectorXd x_aug = VectorXd(n_aug_);
    
    //augmented state covariance
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
    
    //sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_ + 1);
    
    //augmented mean state
    x_aug.fill(0.0);
    x_aug.head(n_x_) = x_;

    //augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    P_aug(5,5) = std_a_*std_a_;
    P_aug(6,6) = std_yawdd_*std_yawdd_;
    
    //square root matrix
    MatrixXd L = P_aug.llt().matrixL();
    
    //augmented sigma points
    Xsig_aug.col(0) = x_aug;
    for(int i = 0; i < n_aug_; i++){
        Xsig_aug.col(i+1) = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
        }
    
    //----------------------------------------------------------------------------
    //   PREDICT SIGMA POINTS
    //----------------------------------------------------------------------------
    
    for (int i = 0; i< 2*n_aug_+1; i++){
        double p_x = Xsig_aug(0,i);
        double p_y = Xsig_aug(1,i);
        double v = Xsig_aug(2,i);
        double yaw = Xsig_aug(3,i);
        double yawd = Xsig_aug(4,i);
        double nu_a = Xsig_aug(5,i);
        double nu_yawdd = Xsig_aug(6,i);
        
        double px_p, py_p;
       
        double sin_yaw = sin(yaw);
        double cos_yaw = cos(yaw);
       
        //verify yaw dot is not a zero so you don't devide by zero
        if(fabs(yawd) > 0.001){
            px_p = p_x + v/yawd * ( sin(yaw + yawd*delta_t) - sin_yaw);
            py_p = p_y + v/yawd * (-cos(yaw + yawd*delta_t) + cos_yaw);
        }
        else {
            px_p = p_x + v*delta_t*cos_yaw;
            py_p = p_y + v*delta_t*sin_yaw;
        }
       
       //add noise to prediction points
        double v_p = v;
        double yaw_p = yaw + yawd*delta_t;
        double yawd_p = yawd;
        
        px_p = px_p + 0.5*nu_a*delta_t*delta_t*cos_yaw;
        py_p = py_p + 0.5*nu_a*delta_t*delta_t*sin_yaw;
        v_p = v_p + delta_t*nu_a;
        yaw_p = yaw_p + 0.5*delta_t*delta_t*nu_yawdd;
        yawd_p = yawd_p + delta_t*nu_yawdd;
        
        Xsig_pred_(0,i) = px_p;
        Xsig_pred_(1,i) = py_p;
        Xsig_pred_(2,i) = v_p;
        Xsig_pred_(3,i) = yaw_p;
        Xsig_pred_(4,i) = yawd_p;
        
        x_ = Xsig_pred_* weights_;
       
       //predicted covarience matrix
        P_.fill(0.0);
        for (int i=0; i<2*n_aug_ + 1;i++){
            VectorXd x_diff = Xsig_pred_.col(i) - x_;
            x_diff(3) = atan2( sin(x_diff(3)), cos(x_diff(3))); // Normalizing the measured value to be within [-pi, pi]
            P_ = P_ + weights_(i)*x_diff*x_diff.transpose();
        }
    }
}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    int n_z = 2;
    MatrixXd Zsig = Xsig_pred_.block(0,0,n_z, 2*n_aug_+1);
    Update(meas_package, Zsig, n_z);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */

void UKF::UpdateRadar(MeasurementPackage meas_package) {   

    int n_z = 3;
    MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_ + 1);
   
    for (int i=0; i < 2*n_aug_+1; i++){
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);
        double v = Xsig_pred_(2,i);
        double yaw = Xsig_pred_(3,i);
        
        double v1 = cos(yaw)*v;
        double v2 = sin(yaw)*v;
        //transform sigma points to measurment space
        Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);
        Zsig(1, i) = atan2(p_y,p_x);
        Zsig(2, i) = (p_x*v1 + p_y*v2)/Zsig(0,i);
    }
    
    Update(meas_package, Zsig, n_z);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement or lidar measurements.
 */

void UKF::Update(MeasurementPackage meas_package, MatrixXd Zsig, int n_z){
   
    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred = Zsig * weights_;
   
    //innovation covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.0);
    for(int i = 0; i < 2*n_aug_+1; i++){
        VectorXd z_diff = Zsig.col(i) - z_pred;
        z_diff(1) = atan2( sin(z_diff(1)), cos(z_diff(1))); // Normalizing the measured value to be within [-pi, pi]
        S = S + weights_(i)*z_diff*z_diff.transpose();
    }
   
    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z,n_z);
    R.fill(0.0);
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
        R << std_radr_*std_radr_, 0, 0,
            0, std_radphi_*std_radphi_, 0,
            0, 0, std_radrd_*std_radrd_;
    }

    if (meas_package.sensor_type_ == MeasurementPackage::LASER){
        R << std_laspx_*std_laspx_, 0,
            0, std_laspy_*std_laspy_;
    }
    
    S = S + R;
   
    //calculate cross correlation matrix
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);
    for (int i = 0; i < 2*n_aug_; i++){
       
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
            z_diff(1) = atan2( sin(z_diff(1)), cos(z_diff(1))); // Normalizing the measured value to be within [-pi, pi]
        }
       
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        x_diff(3) = atan2( sin(x_diff(3)), cos(x_diff(3))); // Normalizing the measured value to be within [-pi, pi]

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    
    VectorXd z = meas_package.raw_measurements_;
    VectorXd z_diff = z - z_pred;
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
        z_diff(1) = atan2( sin(z_diff(1)), cos(z_diff(1))); // Normalizing the measured value to be within [-pi, pi]
    }
   
    //update state mean and covariance matrix
    x_ = x_ + K*z_diff;
    P_ = P_ - K*S*K.transpose();
}
