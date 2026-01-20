#include <Eigen/Core>
#include <iostream> 
#include <fstream>
#include <sstream>
#include <vector>

class MagCal {
private:
    double theta;
    double alpha;
    Eigen::Matrix3d R1;
    Eigen::Matrix3d R2;
    Eigen::Vector3d V;
    Eigen::Vector3d S;
    Eigen::Matrix3d S_diag;
    Eigen::Vector3d axis;
    // test data
    Eigen::Matrix3d R11 = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d R22 = Eigen::Matrix3d::Zero();



public:
    MagCal();
    void calibrateMagnetometer(const Eigen::Matrix3d& R1, const Eigen::Matrix3d& R2, const Eigen::Vector3d& V, const Eigen::Matrix3d S_diag, const Eigen::Vector3d& rawMagData, Eigen::Vector3d& calibratedMagData);
    void calculateCalibrationParameters();
    void update();







};