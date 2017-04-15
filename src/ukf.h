#ifndef UKF_H
#define UKF_H
#include "eigen3/Eigen/Dense"
#include "measurement_package.h"
#include "ground_truth_package.h"
#include <vector>
#include "tools.h"
#include "templates.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::range_error;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* previous measurement
  long previous_timestamp_;
  MeasurementPackage previous_measurement_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* augmented sigma points matrix
  MatrixXd Xsig_aug_;

  ///* time when the state is true, in us
  long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension (includes the two mean noise values)
  int n_aug_;

  ///* Number of sigma points (rule of thumb: 2*n_x_ + 1)
  int n_sigma_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* the current NIS for radar
  double NIS_radar_;

  ///* the current NIS for laser
  double NIS_lidar_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   * @param gt_package The ground truth of the state x at measurement time
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

  void GenerateAugmentedSigmaPoints(const VectorXd &x, const MatrixXd &P, MatrixXd &Xsig_aug);

  /**
   * Predict new sigma points using the CTRV process model.
   *
   * @param delta_t The change in time (in seconds) between the last measurement and this one
   * @param Xsig_aug Reference to the augmented sigma points
   * @param Xsig_pred Reference to the predicted sigma points
   */
  void SigmaPointPrediction(double delta_t, const MatrixXd &Xsig_aug, MatrixXd &Xsig_pred);

  void PredictMeanAndCovariance(const MatrixXd &Xsig_pred, VectorXd &x, MatrixXd &P);

  void PredictLidarMeasurement(const MatrixXd &Xsig_pred, MatrixXd &Ysig, VectorXd &y_pred, MatrixXd &S);

  void PredictRadarMeasurement(const MatrixXd &Xsig_pred, MatrixXd &Zsig, VectorXd &z_pred, MatrixXd &S);

  void UpdateState(const VectorXd &z, const VectorXd &z_pred, const MatrixXd &S,
                   const MatrixXd &Xsig_pred, const MatrixXd &Zsig, VectorXd &x, MatrixXd &P);

  auto InitialiseCovariance() -> MatrixXd;

  auto InitialiseStateVector(MeasurementPackage meas_package) -> VectorXd;

  void DebugPrintState();
};

#endif /* UKF_H */
