#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;


/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;
    previous_timestamp_ = 0;

    // initializing matrices
    R_laser_ = MatrixXd(2, 2);  // laser sensor measurement covariance matrix
    R_radar_ = MatrixXd(3, 3);  // radar sensor measurement covariance matrix
    H_laser_ = MatrixXd(2, 4);  // laser sensor measurement matrix
    Hj_ = MatrixXd(3, 4);       // radar sensor Jacobian matrix

    // defining the constant matrices
    // Laser measurement matrix
    H_laser_ << 1,0,0,0,
          0,1,0,0;
    R_radar_ << 0.09, 0, 0,
          0, 0.0009, 0,
        0, 0, 0.09;
    R_laser_ << 0.0225, 0,
         0, 0.0225;


    double noise_ax = 9;   // x uncertainty added
    double noise_ay = 9;   // y uncertainty added
    ekf_.x_ = VectorXd(4);    // object state
    ekf_.P_ = MatrixXd(4,4);  // object covariance matrix
    // easier to update F matrix from its identity matrix form
    ekf_.F_ = MatrixXd::Identity(4,4);  // state transition matrix
    //ekf_.F_(0,2) = 1;
    //ekf_.F_(1,3) = 1;
    ekf_.Q_ = MatrixXd(4,4);  // process covariance matrix
    //measurement matrix
    ekf_.H_ = MatrixXd(2, 4);
    ekf_.H_ << 1, 0, 0, 0,
            0, 1, 0, 0;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    /*****************************************************************************
     *  Initialization
     ****************************************************************************/

    if (fabs(measurement_pack.raw_measurements_[0]) < 0.0001 or fabs(measurement_pack.raw_measurements_[1]) < 0.0001) {
        cout << "Error in measuring data, skipping it";
        return;
    }
    if (!is_initialized_) {
        /**
        TODO:
          * Initialize the state ekf_.x_ with the first measurement.
           * Create the covariance matrix.
          * Remember: you'll need to convert radar from polar to cartesian coordinates.
        */
        // first measurement
        ekf_.x_ = VectorXd(4);

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            /**
            Convert radar from polar to cartesian coordinates and initialize state.
            */
            double x_cart = measurement_pack.raw_measurements_[0] * cos(measurement_pack.raw_measurements_[1]);
            double y_cart = measurement_pack.raw_measurements_[0] * sin(measurement_pack.raw_measurements_[1]);

            if (x_cart == 0 or y_cart == 0){
                cout << "Error in initializing state matrix";
                return;
            }
            //Initialize the state ekf_.x_ with the first measurement
            ekf_.x_ << x_cart, y_cart, 0, 0;


        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            /**
            Initialize state.
            */
            double x = measurement_pack.raw_measurements_[0];
            double y = measurement_pack.raw_measurements_[1];
            //Initialize the state ekf_.x_ with the first measurement
            ekf_.x_ << x, y, 0, 0;
        }
        // The first measurement is reliable so use the
        // identity matrix as the initial covariance matrix
        ekf_.P_ =  MatrixXd(4,4);
        ekf_.P_ << MatrixXd::Identity(4,4);
        ekf_.P_(3,3)=1000;
        ekf_.P_(2,2)=1000;
        //cout << "EKF initial state matrix is: "  << endl;
        //cout << ekf_.x_  << endl;

        //cout << "EKF initial covariance matrix is: "  << endl;
        //cout << ekf_.P_  << endl;

        // done initializing, no need to predict or update
        is_initialized_ = true;
        previous_timestamp_ = measurement_pack.timestamp_;
        return;
    }

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/


    //cout << measurement_pack.raw_measurements_ << endl;
    //cout << "previous_timestamp_ is:  "<< previous_timestamp_<<endl;
    // divide dt to get seconds units
    double dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;
    //cout << "delta time is:  "<< dt<<endl;
    //cout << "measurement_pack.timestamp_ is:  "<< measurement_pack.timestamp_<<endl;
    // update F matrix to include the delta time, dt
    ekf_.F_(0,2) = dt;
    ekf_.F_(1,3) = dt;

    // Update the process covariance matrix, Q to include time delta


    double dt_2 = dt * dt;
    double dt_3 = dt_2 * dt;
    double
            dt_4 = dt_3 * dt;

    //set the process covariance matrix Q
    ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
            0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
            dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
            0, dt_3/2*noise_ay, 0, dt_2*noise_ay;


    ekf_.Predict();

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    /**
     TODO:
       * Use the sensor type to perform the update step.
       * Update the measurement state and measurement covariance matrices.
     */

    //measurement update
    //cout << " The size of measurement_pack.raw_measurements_ is:  " <<  measurement_pack.raw_measurements_.size() << endl;
    if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        // laser updates
        ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_, ekf_.H_, R_laser_, ekf_.Q_);
        ekf_.Update(measurement_pack.raw_measurements_);
    } else {
        // radar updates
        Hj_ = tools.CalculateJacobian(ekf_.x_);
        ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_, Hj_, R_radar_, ekf_.Q_);
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    }

    // print the output
    //cout << "x_ = " << ekf_.x_ << endl;
    //cout << "P_ = " << ekf_.P_ << endl;
}
