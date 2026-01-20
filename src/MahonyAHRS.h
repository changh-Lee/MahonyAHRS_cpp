//=============================================================================================
// MahonyAHRS.h
//=============================================================================================
//
// Madgwick's implementation of Mayhony's AHRS algorithm.
// See: http://www.x-io.co.uk/open-source-imu-and-ahrs-algorithms/
//
// Date			Author			Notes
// 29/09/2011	SOH Madgwick    Initial release
// 02/10/2011	SOH Madgwick	Optimised for reduced CPU load
//
//=============================================================================================
#ifndef Mahony_AHRS_H
#define Mahony_AHRS_H
#include <math.h>
#include <iostream>
#include <Eigen/Core>

// using namespace Eigen;
// using namespace std;
//---------------------------------------------------------------------------------------------------
// Definitions

constexpr float sampleFreq = 100.0f;				// imu sample frequency (Hz)
constexpr float twoKpDef = (2.0f * 0.5f);			// 2 * proportional gain
constexpr float twoKiDef = (2.0f * 0.0f);			// 2 * integral gain
constexpr float MagneticField[] = {37693.9f, -2256.4f, 25760.2f};	// Earth’s magnetic field(nT) at (113.976274 E,22.594915 N), only need field direction
constexpr int NUM_SAMPLES = 50; 					// a_avg sample amount 11520 100Hz

//--------------------------------------------------------------------------------------------
// Variable declaration

class Mahony {
private:
	float twoKp;		// 2 * proportional gain (Kp)
	float twoKi;		// 2 * integral gain (Ki)
	float q0, q1, q2, q3;	// quaternion of sensor frame relative to auxiliary frame
	float integralFBx, integralFBy, integralFBz;  // integral error terms scaled by Ki
	float invSampleFreq;
	float roll, pitch, yaw;
	char anglesComputed;

	typedef struct {
		float x;
		float y;
		float z;
	} Vector3D;

	typedef struct {
		float b0, b1, b2, a1, a2; // Filter coefficients
		float prev_x[2]; // Previous X-axis input samples
		float prev_y[2]; // Previous Y-axis input samples
		float prev_z[2]; // Previous Z-axis input samples
		float prev_out_x[2]; // Previous X-axis output samples
		float prev_out_y[2]; // Previous Y-axis output samples
		float prev_out_z[2]; // Previous Z-axis output samples
	} SecondOrderFilter3D_t;

	
	void Acc_iirLPF(float ax, float ay, float az, Vector3D *af);
	void SecondOrderFilter3D_init(SecondOrderFilter3D_t* filter, float b0, float b1, float b2, float a1, float a2);
	void SecondOrderFilter3D_update_axis(SecondOrderFilter3D_t* filter, float* input, float* output, char axis);
	void SecondOrderFilter3D_update(SecondOrderFilter3D_t* filter, float ax, float ay, float az, float* out_x, float* out_y, float* out_z);

	float gain1(float ax, float ay, float az);
	float calculate_average_acceleration(Vector3D acceleration_data[], int num_samples);
	float gain2(Vector3D m_perp, Vector3D x, Vector3D af, Vector3D mc);
	float Vector3D_magnitude(Vector3D v);
	float Vector3D_dot_product(Vector3D v, Vector3D w);
	float Vector3D_angle_between(Vector3D v, Vector3D w);


	static float invSqrt(float x);
	void computeAngles();
	Vector3D Acc_buff[NUM_SAMPLES];
//-------------------------------------------------------------------------------------------
// Function declarations

public:
	Mahony();
	void begin(float sampleFrequency) { invSampleFreq = 1.0f / sampleFrequency; }
	void update(float gx, float gy, float gz, float ax, float ay, float az, float mx, float my, float mz);
	void updateIMU(float gx, float gy, float gz, float ax, float ay, float az);
	void Costanzi_update(float gx, float gy, float gz, float ax, float ay, float az, float mx, float my, float mz);
	float getRoll() {
		if (!anglesComputed) computeAngles();
		return roll * 57.29578f;
	}
	float getPitch() {
		if (!anglesComputed) computeAngles();
		return pitch * 57.29578f;
	}
	float getYaw() {
		if (!anglesComputed) computeAngles();
		return yaw * 57.29578f + 180.0f;
	}
	float getRollRadians() {
		if (!anglesComputed) computeAngles();
		return roll;
	}
	float getPitchRadians() {
		if (!anglesComputed) computeAngles();
		return pitch;
	}
	float getYawRadians() {
		if (!anglesComputed) computeAngles();
		return yaw;
	}
};

#endif
