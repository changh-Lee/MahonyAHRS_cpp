#include "magcal.h"
/**
 * @brief This class is used to calibrate magnetometer data.
 * 
*/


/**
 * @brief Default constructor for the MagCal class.
 * 
 * Initializes the member variables theta, axis, alpha, V, S, R1, R2, S_diag, R11, and R22.
 * The member variables theta, axis, alpha, V, and S are assigned with specific values.
 * The member variables R1, R2, S_diag, R11, and R22 are initialized as zero matrices.
 */
MagCal::MagCal()
{

    std::ifstream planefitting_file("data/PlaneFitting_param.csv");
    std::ifstream halir_file("data/Halir_param.csv");

    if (planefitting_file.is_open())
    {
        std::string line;
        std::getline(planefitting_file, line);
        theta = std::stod(line);

        std::getline(planefitting_file, line);
        std::stringstream ss(line);
        for (int i = 0; i < 3; ++i)
        {
            std::string val;
            std::getline(ss, val, ',');
            axis[i] = std::stod(val);
        }
    }

    if (halir_file.is_open())
    {
        std::string line;
        std::getline(halir_file, line);
        alpha = std::stod(line);

        std::getline(halir_file, line);
        std::stringstream ss(line);
        for (int i = 0; i < 3; ++i)
        {
            std::string val;
            std::getline(ss, val, ',');
            V[i] = std::stod(val);
        }

        std::getline(halir_file, line);
        ss = std::stringstream(line);
        for (int i = 0; i < 3; ++i)
        {
            std::string val;
            std::getline(ss, val, ',');
            S[i] = std::stod(val);
        }
    }


    // theta = 0.5755823117859792;// rotation angle in radians
    // axis << 0.25843322, -0.96602912, 0.0;// rotation axis
    // alpha = -0.181879999537509; // rotation angle in radians
    // V << 1.358674116473507, 1.000880842226701e+02, 0.0;//offset
    // S << 0.022075627441017, 0.039882623547207, 1.0;// scale factor

    R1 = Eigen::Matrix3d::Zero();
    R2 = Eigen::Matrix3d::Zero();
    S_diag = Eigen::Matrix3d::Zero();
    R11 = Eigen::Matrix3d::Zero();
    R22 = Eigen::Matrix3d::Zero();

}

/**
 * Calibrates the magnetometer data using the provided calibration parameters.
 * 
 * @param R The rotation matrix representing the sensor orientation.
 * @param V The offset vector representing the sensor bias.
 * @param S The scale vector representing the sensor scale factor.S = [1/a, 1/b, 1/c].'
 * @param rawMagData magdata fitted to the slope.
 * @param calibratedMagData The calibrated magnetometer data (output parameter).
 */
void MagCal::calibrateMagnetometer(const Eigen::Matrix3d& R1, const Eigen::Matrix3d& R2, const Eigen::Vector3d& V, const Eigen::Matrix3d S_diag, const Eigen::Vector3d& rawMagData, Eigen::Vector3d& calibratedMagData)
{
    calibratedMagData = R2 * S_diag * R2.transpose() * (R1 * rawMagData - V);
}

// receive and translate theta, axis, alpha, V, S, into R1, R2, V, S_diag
void MagCal::calculateCalibrationParameters()
{
    double c = cos(theta);
    double s = sin(theta);
    double t = 1 - c;
    double x = axis[0];
    double y = axis[1];
    double z = axis[2];

    R1 << t*x*x + c, t*x*y - s*z, t*x*z + s*y,
          t*x*y + s*z, t*y*y + c, t*y*z - s*x,
          t*x*z - s*y, t*y*z + s*x, t*z*z + c;

    R2 << cos(alpha), -sin(alpha), 0,
          sin(alpha), cos(alpha), 0,
          0, 0, 1;

    S_diag = S.asDiagonal();

    R11 << 0.84963662, -0.04022538, -0.52583223,
          -0.04022538, 0.98923886, -0.14067124,
          0.52583223, 0.14067124, 0.83887548;

    R22 << 0.983505378823939, 0.180878881642882, 0,
          -0.180878881642882, 0.983505378823939, 0,
          0, 0, 1;

    std::cout << "R1:\n" << R1 << "\n";
    std::cout << "R2^T:\n" << R2.transpose() << "\n";
    // std::cout << "R11:\n" << R11 << "\n";
    // std::cout << "R22:\n" << R22 << "\n";
    std::cout << "V:\n" << V << "\n";
    std::cout << "S:\n" << S_diag << "\n";

}

/**
 * Updates the magnetometer calibration and saves the calibrated data to a CSV file.
 * Reads raw magnetometer data from a CSV file, calibrates each data point, and writes the calibrated data to a new CSV file.
 * 
 * @param None
 * @return None
 */
void MagCal::update(){
  // Read raw magnetometer data from CSV file
    std::string filename = "data/m_projection_lse.csv";//F:\\VS_Code\\c++\\MahonyAHRS_c++\\src\\m_projection_lse.csv
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
    }
    std::cout << "Successfully reading raw magnetometer data from " << filename << std::endl;

    std::vector<Eigen::Vector3d> raw_data, calibrated_data;
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string field;
        Eigen::Vector3d data;
        int i = 0;
        while (std::getline(ss, field, ',')) {
        data(i) = std::stof(field);
        i++;
        }
        raw_data.push_back(data);
    }

    file.close();

    // Calibrate each data point
    calibrated_data.reserve(raw_data.size());  // Pre-allocate memory for efficiency
    for (const auto& raw_mag : raw_data) {
        Eigen::Vector3d calibrated_mag;
        this->calibrateMagnetometer(R1, R2, V, S_diag, raw_mag, calibrated_mag);
        // std::cout<<"calibrated magdata:\n"<<calibrated_mag<<"\n";
        calibrated_data.push_back(calibrated_mag);
    }

    // Write calibrated data to a new CSV file
    std::string output_filename = "data/m_calibrated.csv";//F:\\VS_Code\\c++\\MahonyAHRS_c++\\src\\calibrated_data_new.csv
    std::ofstream output_file(output_filename);

    if (!output_file.is_open()) {
    std::cerr << "Error: Could not open file " << output_filename << std::endl;
    }

    for (const auto& calibrated_mag : calibrated_data) {
    output_file << calibrated_mag(0) << "," << calibrated_mag(1) << "," << calibrated_mag(2) << "\n";
    }

    output_file.close();

    std::cout << "Calibrated data saved to " << output_filename << std::endl;

}


int main(){
    MagCal MagCal;
    MagCal.calculateCalibrationParameters();
    MagCal.update();

    return 0;
}
