#include "kalman_filter.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"
using namespace Eigen;
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

KalmanFilter::KalmanFilter() {}
KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(const VectorXd &x_in, const MatrixXd &P_in, const MatrixXd &F_in,
                        const MatrixXd &H_in, const MatrixXd &R_in, const MatrixXd &Q_in) {
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    H_ = H_in;
    R_ = R_in;
    Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */

        x_ = F_ * x_;
        //std::cout << " The predict x state is "  <<  x_<< endl;
        MatrixXd Ft = F_.transpose();
        P_ = F_ * P_ * Ft + Q_;
        //std::cout << " The covariance matrix, P, is "  <<  P_<< endl;
}
// update(MeasurementModel&, VectorXd&)
void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
    // z is the measured value of the sensors
    H_ = MatrixXd(2, 4);
    H_ << 1, 0, 0, 0,
            0, 1, 0, 0;
    //std::cout << "x_ is: " << x_ << std::endl;
    //std::cout << "H_ in laser is: " << H_ << std::endl;

    VectorXd z_pred = H_ * x_;
    VectorXd y = z - z_pred;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_ ;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);

    long long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
    //std::cout << "P_ new estimate is: " << P_ << std::endl;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

  //TODO:
   //  * update the state by using Extended Kalman Filter equations

    // z is the measured value of the sensors

    /**
            Convert radar from polar to cartesian coordinates and initialize state.
            */

    //double rho = sqrt((x_[0]*x_[0] + x_[1]*x_[1]));
    //double phi = atan(x_[1]/x_[0]);
    //double phi = atan2(x_[1],x_[0]);
    //double rho_dot= (x_[0]*x_[2] + x_[1]*x_[3])/sqrt((x_[0]*x_[0] + x_[1]*x_[1]));

    double rho = sqrt((x_(0)*x_(0) + x_(1)*x_(1)));
    //double phi = atan(x_[1]/x_[0]);
    double phi = atan2(x_(1),x_(0));
    double rho_dot= (x_(0)*x_(2) + x_(1)*x_[3])/sqrt((x_(0)*x_(0) + x_(1)*x_(1)));


    //out << " rho, phi, rho_dot  " << rho << phi << rho_dot<< ::endl;
    MatrixXd  z_pred_radar(3, 1);
    z_pred_radar << rho, phi, rho_dot;
    VectorXd y = z - z_pred_radar;


    MatrixXd Ht = H_.transpose();

    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    //std::cout << "x_ laser new estimate is: " << std::endl;
    //std::cout << x_<< std::endl;
    size_t x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}

