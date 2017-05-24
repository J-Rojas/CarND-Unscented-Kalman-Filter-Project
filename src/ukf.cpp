#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

VectorXd polar_to_cartesian(const VectorXd& in) {
  VectorXd cartesian(4);

  double phi = in(1);
  double cosangle = cos(phi);
  double sinangle = sin(phi);
  double range = in(0);
  double range_rate = in(2);
  double x = cosangle * range;
  double y = sinangle * range;
  double rate_x = cosangle * range_rate;
  double rate_y = sinangle * range_rate;

  cartesian << x, y, rate_x, rate_y;

  return cartesian;
}

VectorXd GenerateWeights(double lambda, int n_size) {
	//set vector for weights
	VectorXd weights(2*n_size+1);
	double weight_0 = lambda/(lambda+n_size);
	weights(0) = weight_0;
	for (int i=1; i<2*n_size+1; i++) {
		double weight = 0.5/(n_size+lambda);
		weights(i) = weight;
	}
	return weights;
}

double NormalizeAngle(double value)
{
	if (fabs(value) > M_PI)
	{
		if (fabs(value) > 100000) {
			std::cout << "Something is wrong " << std::endl;
		}
		value -= round(value / (2.0 * M_PI)) * (2.0 * M_PI);
	}
	return value;
}

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
	VectorXd x(5);
  x_ = x;

  // initial covariance matrix
	MatrixXd P(5, 5);
	P_ = P;
	P_ << 1, 0, 0, 0, 0,
				0, 1, 0, 0, 0,
				0, 0, 1, 0, 0,
				0, 0, 0, 1, 0,
				0, 0, 0, 0, 1;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.55;

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

	// Dimension of state vector
	n_x_ = 5;

	// Dimension of augmented state vector
	n_aug_ = 7;

	// Lambda design parameter
	lambda_ = 3 - n_aug_;

	weights_ = GenerateWeights(lambda_, n_aug_);

	is_initialized_ = false;

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

  if (!is_initialized_) {

    //initialize state and prediction
    switch (meas_package.sensor_type_) {
      case MeasurementPackage::SensorType::LASER:
          x_ << meas_package.raw_measurements_(0),
                meas_package.raw_measurements_(1),
                0,
                0,
                0;

          is_initialized_ = true;
        break;
      case MeasurementPackage::SensorType::RADAR:
          //convert from polar coordinates
          VectorXd cart = polar_to_cartesian(meas_package.raw_measurements_);

          x_ << cart(0),
                cart(1),
                cart(2),
                cart(3),
                0;

          is_initialized_ = true;
        break;
    }

  } else {

    double delta = (meas_package.timestamp_ - time_us_) / 1000000.0;

	  Prediction(delta);

    switch (meas_package.sensor_type_) {
      case MeasurementPackage::SensorType::LASER:
        if (use_laser_)
          UpdateLidar(meas_package);
        break;
      case MeasurementPackage::SensorType::RADAR:
        if (use_radar_)
          UpdateRadar(meas_package);
        break;
    }

  }

  time_us_ = meas_package.timestamp_;

}

MatrixXd GenerateStateSigmaPoints(UKF& ukf) {

  //create sigma point matrix
  MatrixXd Xsig(ukf.n_aug_, 2 * ukf.n_aug_ + 1);

  //create augmented mean vector
  VectorXd x_aug(ukf.n_aug_);

  //create augmented state covariance
  MatrixXd P_aug(ukf.n_aug_, ukf.n_aug_);

  //calculate square root of P
  MatrixXd A = ukf.P_.llt().matrixL();

  //create augmented mean state
  x_aug.head(5) = ukf.x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = ukf.P_;
  P_aug(5,5) = ukf.std_a_ * ukf.std_a_;
  P_aug(6,6) = ukf.std_yawdd_ * ukf.std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig.col(0)  = x_aug;
  for (int i = 0; i< ukf.n_aug_; i++)
  {
    Xsig.col(i+1)       = x_aug + sqrt(ukf.lambda_+ukf.n_aug_) * L.col(i);
    Xsig.col(i+1+ukf.n_aug_) = x_aug - sqrt(ukf.lambda_+ukf.n_aug_) * L.col(i);
  }

  return Xsig;

}

MatrixXd GenerateMeasurementSigmaPoints(UKF& ukf) {

	int n_z = 3;

	//create matrix for sigma points in measurement space
	MatrixXd Zsig(n_z, 2 * ukf.n_aug_ + 1);

	//transform sigma points into measurement space
	for (int i = 0; i < 2 * ukf.n_aug_ + 1; i++) {  //2n+1 simga points

		// extract values for better readibility
		double p_x = ukf.Xsig_pred_(0,i);
		double p_y = ukf.Xsig_pred_(1,i);
		double v  = ukf.Xsig_pred_(2,i);
		double yaw = ukf.Xsig_pred_(3,i);

		double v1 = cos(yaw)*v;
		double v2 = sin(yaw)*v;

		// measurement model
		double norm = sqrt(p_x*p_x + p_y*p_y);

		Zsig(0,i) = norm;                        //r
		Zsig(1,i) = atan2(p_y,p_x);              //phi
		Zsig(2,i) = norm > 0.0001 ? (p_x*v1 + p_y*v2 ) / norm : 0.0;   //r_dot
	}

	return Zsig;

}


MatrixXd UKF::TransformedSigmaPoints(MatrixXd& Xsig_aug, double delta_t) {

  MatrixXd Xsig_pred(n_x_, 2 * n_aug_ + 1);

  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x      = Xsig_aug(0,i);
    double p_y      = Xsig_aug(1,i);
    double v        = Xsig_aug(2,i);
    double yaw      = Xsig_aug(3,i);
    double yawd     = Xsig_aug(4,i);
    double nu_a     = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;
	  double v_p = v;
	  double yaw_p = yaw + yawd * delta_t;
	  double yawd_p = yawd;

	  //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v/yawd * ( sin(yaw_p) - sin(yaw) );
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw_p) );
    }
    else {
      px_p = p_x + v * delta_t * cos(yaw);
      py_p = p_y + v * delta_t * sin(yaw);
    }

    //add noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
	  yawd_p = yawd_p + nu_yawdd * delta_t;

	  //write predicted sigma point into correct column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }

  return Xsig_pred;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  MatrixXd sigma = GenerateStateSigmaPoints(*this);

#ifdef PRINT
	std::cout << "sigma " << std::endl << sigma << std::endl;
#endif

  Xsig_pred_ = TransformedSigmaPoints(sigma, delta_t);

#ifdef PRINT
	std::cout << "tran(sigma)" << std::endl << Xsig_pred_ << std::endl;
#endif

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

#ifdef PRINT
	std::cout << "Predicted mean" << std::endl << x_ << std::endl;
#endif

	//predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
	  x_diff(3) = NormalizeAngle(x_diff(3));

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }

#ifdef PRINT
	std::cout << "Predicted covariance" << std::endl << P_ << std::endl;
#endif

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.

  LIDAR measurement is in linear because it is in cartesian space, unlike radar which is in
   polar (non-linear) coordinates. For LIDAR measurements, we use the standard Kalman filter
   equations.

  */

	VectorXd z = meas_package.raw_measurements_;
	VectorXd x = x_;
	MatrixXd H(2, n_x_);
	MatrixXd R(2, 2);
	H << 1, 0, 0, 0, 0,
			 0, 1, 0, 0, 0;
	R << std_laspx_ * std_laspx_, 0,
			 0, std_laspy_ * std_laspy_;

	VectorXd z_pred = H * x;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H.transpose();
	MatrixXd S = H * P_ * Ht + R;
	MatrixXd K = P_ * Ht * S.inverse();

	x_ = x_ + K * y;

	x_(3) = NormalizeAngle(x_(3));

	P_ = (Eigen::Matrix<double, 5, 5>::Identity() - K * H) * P_;

	// the apriori estimate of the position is used
	NIS_laser_ = Tools::CalculateNIS(y, S);

}

void UKF::PredictRadarMeasurement(MatrixXd& Zsig, VectorXd& z_pred, MatrixXd& S_out) {

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

	//mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
	  z_diff(1) = NormalizeAngle(z_diff(1));

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0, std_radrd_*std_radrd_;
  S_out = S + R;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

	/**
  Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

	//calculate cross correlation matrix
	int n_z = 3;
	VectorXd z = meas_package.raw_measurements_;

	MatrixXd S(n_z,n_z);
	VectorXd z_pred(n_z);

	//Generate Sigma Points in Measurement Space
	MatrixXd Zsig = GenerateMeasurementSigmaPoints(*this);

	//Calculate predicted sigma points in measurement space
	PredictRadarMeasurement(Zsig, z_pred, S);

	MatrixXd Tc(n_x_, n_z);
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		VectorXd xv = Xsig_pred_.col(i) - x_;

		xv(3) = NormalizeAngle(xv(3));

		VectorXd zv = Zsig.col(i) - z_pred;

		zv(1) = NormalizeAngle(zv(1));

		Tc += weights_(i) * xv * zv.transpose();
	}

	//calculate Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	//update state mean and covariance matrix
	VectorXd z_diff = (z - z_pred);

	z_diff(1) = NormalizeAngle(z_diff(1));

	x_ = x_ + K * z_diff;
	P_ = P_ - K * S * K.transpose();

	// update nis metric
	NIS_radar_ = Tools::CalculateNIS(z_diff, S);

}
