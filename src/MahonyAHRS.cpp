//=============================================================================================
// MahonyAHRS.c
//=============================================================================================
//
// Madgwick's implementation of Mayhony's AHRS algorithm.
// See: http://www.x-io.co.uk/open-source-imu-and-ahrs-algorithms/
//
// From the x-io website "Open-source resources available on this website are
// provided under the GNU General Public Licence unless an alternative licence
// is provided in source."
//
// Date			Author			Notes
// 29/09/2011	SOH Madgwick    Initial release
// 02/10/2011	SOH Madgwick	Optimised for reduced CPU load
//
// Algorithm paper:
// http://ieeexplore.ieee.org/xpl/login.jsp?tp=&arnumber=4608934&url=http%3A%2F%2Fieeexplore.ieee.org%2Fstamp%2Fstamp.jsp%3Ftp%3D%26arnumber%3D4608934
//
//=============================================================================================

//-------------------------------------------------------------------------------------------
// Header files

#include "MahonyAHRS.h"
#include "magcal.h"

//-------------------------------------------------------------------------------------------
// Definitions

// #define DEFAULT_SAMPLE_FREQ	512.0f	// sample frequency in Hz
// #define twoKpDef	(2.0f * 0.5f)	// 2 * proportional gain
// #define twoKiDef	(2.0f * 0.0f)	// 2 * integral gain


//============================================================================================
// Functions

//-------------------------------------------------------------------------------------------
// AHRS algorithm update

Mahony::Mahony()
{
	twoKp = twoKpDef;	// 2 * proportional gain (Kp)
	twoKi = twoKiDef;	// 2 * integral gain (Ki)
	q0 = 1.0f;
	q1 = 0.0f;
	q2 = 0.0f;
	q3 = 0.0f;
	integralFBx = 0.0f;
	integralFBy = 0.0f;
	integralFBz = 0.0f;
	anglesComputed = 0;
	invSampleFreq = 1.0f / sampleFreq;
}

/**
 * Updates the Mahony AHRS algorithm with new sensor data.
 *
 * @param gx The gyroscope X-axis reading. rad/s
 * @param gy The gyroscope Y-axis reading.
 * @param gz The gyroscope Z-axis reading.
 * @param ax The accelerometer X-axis reading.
 * @param ay The accelerometer Y-axis reading.
 * @param az The accelerometer Z-axis reading.
 * @param mx The magnetometer X-axis reading.
 * @param my The magnetometer Y-axis reading.
 * @param mz The magnetometer Z-axis reading.
 */
void Mahony::Costanzi_update(float gx, float gy, float gz, float ax, float ay, float az, float mx, float my, float mz)
{
	float recipNorm;
	float q0q0, q0q1, q0q2, q0q3, q1q1, q1q2, q1q3, q2q2, q2q3, q3q3;
	float _m_perp;
	float halfvx, halfvy, halfvz, halfwx, halfwy, halfwz;
	float halfex, halfey, halfez;
	float qa, qb, qc;
	float k1, k2;
	Vector3D x = {0.0, 0.0, 0.0};// magnetic north(fixed-frame x-coordinate) in body frame
	Vector3D af = {0.0, 0.0, 0.0};// filtered acc
	Vector3D mc = {0.0, 0.0, 0.0};// calibrated magnetic field
	Vector3D m_perp = {0.0, 0.0, 0.0};// projection of m^c orthogonal to af

	// Use IMU algorithm if magnetometer measurement invalid(avoids NaN in magnetometer normalisation)
	if((mx == 0.0f) && (my == 0.0f) && (mz == 0.0f)) {
		updateIMU(gx, gy, gz, ax, ay, az);
		return;
	}

	// Compute feedback only if accelerometer measurement valid(avoids NaN in accelerometer normalisation)
	if(!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f))) {
		// gain k1
		k1 = gain1(ax, ay, az);

		// Normalise accelerometer measurement
		recipNorm = invSqrt(ax * ax + ay * ay + az * az);
		ax *= recipNorm;
		ay *= recipNorm;
		az *= recipNorm;

		// Normalise filtered acc
		Acc_iirLPF(ax, ay, az, &af);// add a second-order-filter to acc
		recipNorm = invSqrt(af.x * af.x + af.y * af.y + af.z * af.z);
		af.x *= recipNorm;
		af.y *= recipNorm;
		af.z *= recipNorm; 

		// Normalise magnetometer measurement
		MagCal MagCal;
		MagCal.calculateCalibrationParameters();
		MagCal.update();
		recipNorm = invSqrt(mx * mx + my * my + mz * mz);//这里需要对磁场进行额外的初始校正
		mx *= recipNorm;
		my *= recipNorm;
		mz *= recipNorm;
		mc.x = mx;
		mc.y = my;
		mc.z = mz; 

		// Auxiliary variables to avoid repeated arithmetic
		q0q0 = q0 * q0;
		q0q1 = q0 * q1;
		q0q2 = q0 * q2;
		q0q3 = q0 * q3;
		q1q1 = q1 * q1;
		q1q2 = q1 * q2;
		q1q3 = q1 * q3;
		q2q2 = q2 * q2;
		q2q3 = q2 * q3;
		q3q3 = q3 * q3;

		// choose magnetic north as a known direction m_perp^c
		_m_perp = af.x*mx +af.y*my + af.z*mz;// norm of m_perp^c
		m_perp.x = mx - af.x*_m_perp;
		m_perp.y = my - af.y*_m_perp;
		m_perp.z = mz - af.z*_m_perp;

		// Estimated direction of gravity and magnetic field
		halfvx = q1q3 - q0q2;
		halfvy = q0q1 + q2q3;
		halfvz = q0q0 - 0.5f + q3q3;
		halfwx = 0.5f - q2q2 - q3q3;// magnetic north => R*x^N = R[1 0 0]^T
        halfwy = q1q2 - q0q3;
        halfwz = q0q2 + q1q3;  
		x.x = halfwx;//magnetic north in body frame
		x.y = halfwy;
		x.z = halfwz;

		//gain k2
		k2 = gain2(m_perp, x, af, mc);

		// Error is sum of cross product between estimated direction and measured direction of field vectors，对应ω_mes = k1(a_f x Rz^N) + k2(m_perp x Rx^N) = error_grav + error_magn
		halfex = k1*(af.y * halfvz - af.z * halfvy) + k2*(m_perp.y * halfwz - m_perp.z * halfwy);//这里加上k1，k2及其定义
		halfey = k1*(af.z * halfvx - af.x * halfvz) + k2*(m_perp.z * halfwx - m_perp.x * halfwz);
		halfez = k1*(af.x * halfvy - af.y * halfvx) + k2*(m_perp.x * halfwy - m_perp.y * halfwx);

		// Compute and apply integral feedback if enabled
		if(twoKi > 0.0f) {
			// integral error scaled by Ki
			integralFBx += twoKi * halfex * invSampleFreq;
			integralFBy += twoKi * halfey * invSampleFreq;
			integralFBz += twoKi * halfez * invSampleFreq;
			gx += integralFBx;	// apply integral feedback
			gy += integralFBy;
			gz += integralFBz;
		} else {
			integralFBx = 0.0f;	// prevent integral windup
			integralFBy = 0.0f;
			integralFBz = 0.0f;
		}

		// Apply proportional feedback
		gx += twoKp * halfex;
		gy += twoKp * halfey;
		gz += twoKp * halfez;
	}

	// Integrate rate of change of quaternion
	gx *= (0.5f * invSampleFreq);		// pre-multiply common factors
	gy *= (0.5f * invSampleFreq);
	gz *= (0.5f * invSampleFreq);
	qa = q0;
	qb = q1;
	qc = q2;
	q0 += (-qb * gx - qc * gy - q3 * gz);
	q1 += (qa * gx + qc * gz - q3 * gy);
	q2 += (qa * gy - qb * gz + q3 * gx);
	q3 += (qa * gz + qb * gy - qc * gx);

	// Normalise quaternion
	recipNorm = invSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
	q0 *= recipNorm;
	q1 *= recipNorm;
	q2 *= recipNorm;
	q3 *= recipNorm;
	anglesComputed = 0;
}

void Mahony::update(float gx, float gy, float gz, float ax, float ay, float az, float mx, float my, float mz)
{
	float recipNorm;
	float q0q0, q0q1, q0q2, q0q3, q1q1, q1q2, q1q3, q2q2, q2q3, q3q3;
	float hx, hy, bx, bz;
	float halfvx, halfvy, halfvz, halfwx, halfwy, halfwz;
	float halfex, halfey, halfez;
	float qa, qb, qc;

	// Use IMU algorithm if magnetometer measurement invalid
	// (avoids NaN in magnetometer normalisation)
	if((mx == 0.0f) && (my == 0.0f) && (mz == 0.0f)) {
		updateIMU(gx, gy, gz, ax, ay, az);
		return;
	}

	// Convert gyroscope degrees/sec to radians/sec
	gx *= 0.0174533f;
	gy *= 0.0174533f;
	gz *= 0.0174533f;

	// Compute feedback only if accelerometer measurement valid
	// (avoids NaN in accelerometer normalisation)
	if(!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f))) {

		// Normalise accelerometer measurement
		recipNorm = invSqrt(ax * ax + ay * ay + az * az);
		ax *= recipNorm;
		ay *= recipNorm;
		az *= recipNorm;

		// Normalise magnetometer measurement
		recipNorm = invSqrt(mx * mx + my * my + mz * mz);
		mx *= recipNorm;
		my *= recipNorm;
		mz *= recipNorm;

		// Auxiliary variables to avoid repeated arithmetic
		q0q0 = q0 * q0;
		q0q1 = q0 * q1;
		q0q2 = q0 * q2;
		q0q3 = q0 * q3;
		q1q1 = q1 * q1;
		q1q2 = q1 * q2;
		q1q3 = q1 * q3;
		q2q2 = q2 * q2;
		q2q3 = q2 * q3;
		q3q3 = q3 * q3;

		// Reference direction of Earth's magnetic field
		hx = 2.0f * (mx * (0.5f - q2q2 - q3q3) + my * (q1q2 - q0q3) + mz * (q1q3 + q0q2));
		hy = 2.0f * (mx * (q1q2 + q0q3) + my * (0.5f - q1q1 - q3q3) + mz * (q2q3 - q0q1));
		bx = sqrtf(hx * hx + hy * hy);
		bz = 2.0f * (mx * (q1q3 - q0q2) + my * (q2q3 + q0q1) + mz * (0.5f - q1q1 - q2q2));

		// Estimated direction of gravity and magnetic field
		halfvx = q1q3 - q0q2;
		halfvy = q0q1 + q2q3;
		halfvz = q0q0 - 0.5f + q3q3;
		halfwx = bx * (0.5f - q2q2 - q3q3) + bz * (q1q3 - q0q2);
		halfwy = bx * (q1q2 - q0q3) + bz * (q0q1 + q2q3);
		halfwz = bx * (q0q2 + q1q3) + bz * (0.5f - q1q1 - q2q2);

		// Error is sum of cross product between estimated direction
		// and measured direction of field vectors
		halfex = (ay * halfvz - az * halfvy) + (my * halfwz - mz * halfwy);
		halfey = (az * halfvx - ax * halfvz) + (mz * halfwx - mx * halfwz);
		halfez = (ax * halfvy - ay * halfvx) + (mx * halfwy - my * halfwx);

		// Compute and apply integral feedback if enabled
		if(twoKi > 0.0f) {
			// integral error scaled by Ki
			integralFBx += twoKi * halfex * invSampleFreq;
			integralFBy += twoKi * halfey * invSampleFreq;
			integralFBz += twoKi * halfez * invSampleFreq;
			gx += integralFBx;	// apply integral feedback
			gy += integralFBy;
			gz += integralFBz;
		} else {
			integralFBx = 0.0f;	// prevent integral windup
			integralFBy = 0.0f;
			integralFBz = 0.0f;
		}

		// Apply proportional feedback
		gx += twoKp * halfex;
		gy += twoKp * halfey;
		gz += twoKp * halfez;
	}

	// Integrate rate of change of quaternion
	gx *= (0.5f * invSampleFreq);		// pre-multiply common factors
	gy *= (0.5f * invSampleFreq);
	gz *= (0.5f * invSampleFreq);
	qa = q0;
	qb = q1;
	qc = q2;
	q0 += (-qb * gx - qc * gy - q3 * gz);
	q1 += (qa * gx + qc * gz - q3 * gy);
	q2 += (qa * gy - qb * gz + q3 * gx);
	q3 += (qa * gz + qb * gy - qc * gx);

	// Normalise quaternion
	recipNorm = invSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
	q0 *= recipNorm;
	q1 *= recipNorm;
	q2 *= recipNorm;
	q3 *= recipNorm;
	anglesComputed = 0;
}

//-------------------------------------------------------------------------------------------
// IMU algorithm update

void Mahony::updateIMU(float gx, float gy, float gz, float ax, float ay, float az)
{
	float recipNorm;
	float halfvx, halfvy, halfvz;
	float halfex, halfey, halfez;
	float qa, qb, qc;

	// Convert gyroscope degrees/sec to radians/sec
	gx *= 0.0174533f;
	gy *= 0.0174533f;
	gz *= 0.0174533f;

	// Compute feedback only if accelerometer measurement valid
	// (avoids NaN in accelerometer normalisation)
	if(!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f))) {

		// Normalise accelerometer measurement
		recipNorm = invSqrt(ax * ax + ay * ay + az * az);
		ax *= recipNorm;
		ay *= recipNorm;
		az *= recipNorm;

		// Estimated direction of gravity
		halfvx = q1 * q3 - q0 * q2;
		halfvy = q0 * q1 + q2 * q3;
		halfvz = q0 * q0 - 0.5f + q3 * q3;

		// Error is sum of cross product between estimated
		// and measured direction of gravity
		halfex = (ay * halfvz - az * halfvy);
		halfey = (az * halfvx - ax * halfvz);
		halfez = (ax * halfvy - ay * halfvx);

		// Compute and apply integral feedback if enabled
		if(twoKi > 0.0f) {
			// integral error scaled by Ki
			integralFBx += twoKi * halfex * invSampleFreq;
			integralFBy += twoKi * halfey * invSampleFreq;
			integralFBz += twoKi * halfez * invSampleFreq;
			gx += integralFBx;	// apply integral feedback
			gy += integralFBy;
			gz += integralFBz;
		} else {
			integralFBx = 0.0f;	// prevent integral windup
			integralFBy = 0.0f;
			integralFBz = 0.0f;
		}

		// Apply proportional feedback
		gx += twoKp * halfex;
		gy += twoKp * halfey;
		gz += twoKp * halfez;
	}

	// Integrate rate of change of quaternion
	gx *= (0.5f * invSampleFreq);		// pre-multiply common factors
	gy *= (0.5f * invSampleFreq);
	gz *= (0.5f * invSampleFreq);
	qa = q0;
	qb = q1;
	qc = q2;
	q0 += (-qb * gx - qc * gy - q3 * gz);
	q1 += (qa * gx + qc * gz - q3 * gy);
	q2 += (qa * gy - qb * gz + q3 * gx);
	q3 += (qa * gz + qb * gy - qc * gx);

	// Normalise quaternion
	recipNorm = invSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
	q0 *= recipNorm;
	q1 *= recipNorm;
	q2 *= recipNorm;
	q3 *= recipNorm;
	anglesComputed = 0;
}

//-------------------------------------------------------------------------------------------
// Fast inverse square-root
// See: http://en.wikipedia.org/wiki/Fast_inverse_square_root

float Mahony::invSqrt(float x)
{
	float halfx = 0.5f * x;
	union { float f; long l; } i;
	i.f = x;
	i.l = 0x5f3759df - (i.l >> 1);
	float y = i.f;
	y = y * (1.5f - (halfx * y * y));
	y = y * (1.5f - (halfx * y * y));
	return y;
}

//-------------------------------------------------------------------------------------------

float Mahony::gain1(float ax, float ay, float az)
{
	float k1_init = 0.2;// acc more reliable than compass k1_init>k2_init
	float k1 = 0.0;
	float Da = 0.0;//absolute percentage error
	float D_threshold = 0.3;
	float D_max = 0.5;
	float delta_D = 0.0, a_avg = 9.81;

	a_avg = calculate_average_acceleration(Acc_buff, NUM_SAMPLES);
	Da = abs((sqrt(ax*ax+ay*ay+az*az) - a_avg) / a_avg);// unnormalized and unfiltered acc
	if (Da <= D_threshold)
		k1 = k1_init;
	else if (Da > D_threshold && Da < D_max){
		delta_D = (D_max - D_threshold);
		k1 = (k1_init/delta_D) * (D_max - Da);
	}
	else
		k1 = 0;

	return k1;
}

float Mahony::gain2(Vector3D m_perp, Vector3D x, Vector3D af, Vector3D mc){
	float k2_init = 0.1;// 参考论文数值
	float k2=0.0;
	float alpha1=0.0, alpha2=0.0, alpha1_th=5.0*0.0174533, alpha2_th=5.0*0.0174533;//rad
	int cont1 = 0, cont2 = 0;
	int cont1_max = 5;
	int cont2_max = 10;
	Vector3D z = {0.0, 0.0, 1.0};
	Vector3D H = {MagneticField[0], MagneticField[1], MagneticField[2]};

	alpha1 = Vector3D_angle_between(x, m_perp);
	alpha2 = abs(Vector3D_angle_between(mc, af) - Vector3D_angle_between(H, z));

	if (alpha1 > alpha1_th || alpha2 > alpha2_th){
		// k2>=0
		if (cont1 < cont1_max){
			k2 = k2_init*(1 - cont1 / cont1_max);
			cont1++;
		}
		else
			k2 = 0;

		if (cont2 != 0) 
			cont2 = 0;
	}
	else{
		// k2<=k2_init
		if (k2 < k2_init){
			k2 += (k2_init - k2)*(cont2/cont2_max);
			cont2++;
		}
		k2 = k2 > k2_init ? k2_init : k2;

		if (cont1 != 0)
			cont1 = 0;
	}

	return k2;
}

// 计算两个向量之间的夹角（弧度）
float Mahony::Vector3D_angle_between(Vector3D v, Vector3D w){
    float v_magnitude = Vector3D_magnitude(v);
    float w_magnitude = Vector3D_magnitude(w);
    float dot = Vector3D_dot_product(v, w);
    return acos(dot / (v_magnitude * w_magnitude));
}

// 计算加速度在一段时间内的平均值
float Mahony::calculate_average_acceleration(Vector3D acceleration_data[], int num_samples)
{
    int i;
	
	float sum_x = 0.0f;
    float sum_y = 0.0f;
    float sum_z = 0.0f;
    
    for (i = 0; i < num_samples; ++i) {
        sum_x += acceleration_data[i].x;
        sum_y += acceleration_data[i].y;
        sum_z += acceleration_data[i].z;
    }

	sum_x /= num_samples;
	sum_y /= num_samples;
    sum_z /= num_samples;

    return sqrtf(sum_x*sum_x + sum_y*sum_y + sum_z*sum_z);
}

// second-order filter
void Mahony::Acc_iirLPF(float ax, float ay, float az, Vector3D *af){

	// IIR Filter Coefficients
	float b0 = 0.067455273889072f;// fc = 10 Hz
	float b1 = 0.134910547778144f;
	float b2 = 0.067455273889072f;
	float a1 = -1.142980502539901f;// a0 = 1
	float a2 = 0.412801598096189f;

	// Initialize the filter
	SecondOrderFilter3D_t filter;
	SecondOrderFilter3D_init(&filter, b0, b1, b2, a1, a2);
	SecondOrderFilter3D_update(&filter, ax, ay, az, &(af->x), &(af->y), &(af->z));

	// printf("Filtered X-axis acceleration: %f\n", out_x);
	// printf("Filtered Y-axis acceleration: %f\n", out_y);
	// printf("Filtered Z-axis acceleration: %f\n", out_z);
}

void Mahony::SecondOrderFilter3D_init(SecondOrderFilter3D_t* filter, float b0, float b1, float b2, float a1, float a2) {
	int i;
	
	filter->b0 = b0;
	filter->b1 = b1;
	filter->b2 = b2;
	filter->a1 = a1;
	filter->a2 = a2;
	
	for (i = 0; i < 2; i++) {
		filter->prev_x[i] = 0.0f;
		filter->prev_y[i] = 0.0f;
		filter->prev_z[i] = 0.0f;
		filter->prev_out_x[i] = 0.0f;
		filter->prev_out_y[i] = 0.0f;
		filter->prev_out_z[i] = 0.0f;
	}
}

void Mahony::SecondOrderFilter3D_update(SecondOrderFilter3D_t* filter, float ax, float ay, float az, float* out_x, float* out_y, float* out_z) {
	SecondOrderFilter3D_update_axis(filter, &ax, out_x, 'X');
	SecondOrderFilter3D_update_axis(filter, &ay, out_y, 'Y');
	SecondOrderFilter3D_update_axis(filter, &az, out_z, 'Z');
}

void Mahony::SecondOrderFilter3D_update_axis(SecondOrderFilter3D_t* filter, float* input, float* output, char axis) {
    float out_tmp;
    float *prev_input= NULL, *prev_output= NULL;

    switch(axis) {
        case 'X':
            prev_input = filter->prev_x;
            prev_output = filter->prev_out_x;
            break;
        case 'Y':
            prev_input = filter->prev_y;
            prev_output = filter->prev_out_y;
            break;
        case 'Z':
            prev_input = filter->prev_z;
            prev_output = filter->prev_out_z;
            break;
        default:
			std::cout<<"Invalid axis: "<<axis<<std::endl;
            return;
    }

    out_tmp = filter->b0 * *input + filter->b1 * prev_input[0] + filter->b2 * prev_input[1] 
              					  - filter->a1 * prev_output[0] - filter->a2 * prev_output[1];//a1*y_k = b1*x_k + b2*x_k-1 + b3*x_k-2 - a2*y_k-1 - a3*y_k-2

    prev_input[1] = prev_input[0];
    prev_input[0] = *input;
    prev_output[1] = prev_output[0];
    prev_output[0] = out_tmp;

    *output = out_tmp;
}

void Mahony::computeAngles()
{
	roll = atan2f(q0*q1 + q2*q3, 0.5f - q1*q1 - q2*q2);
	pitch = asinf(-2.0f * (q1*q3 - q0*q2));
	yaw = atan2f(q1*q2 + q0*q3, 0.5f - q2*q2 - q3*q3);
	anglesComputed = 1;
}

// 计算向量的模
float Mahony::Vector3D_magnitude(Vector3D v){
    return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
}

// 计算向量的点积
float Mahony::Vector3D_dot_product(Vector3D v, Vector3D w){
    return v.x * w.x + v.y * w.y + v.z * w.z;
}

//============================================================================================
// END OF CODE
//============================================================================================
