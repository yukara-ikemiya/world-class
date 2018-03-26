//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise,
// Copyright 2018 Yukara Ikemiya
//
// Band-aperiodicity estimation on the basis of the idea of D4C.
//-----------------------------------------------------------------------------

#include "d4c.hpp"

#include <cmath>
#include <numeric>
#include <algorithm>
#include <omp.h>

#include "world_common.hpp"
#include "world_constantnumbers.hpp"
#include "world_matlabfunctions.hpp"

namespace
{
using std::copy;
using std::generate;
using std::accumulate;
using std::for_each;
}


namespace world_class
{

D4COption::D4COption()
	: threshold(world::kThreshold)
{}

D4C::D4C(int fs)
{
	prepareForD4c(fs);
}

D4C::D4C(int fs, const D4COption &option)
{
	option_.threshold = option.threshold;
	prepareForD4c(fs);
}

D4C::~D4C()
{
	for (int i = 0; i < num_thread_; ++i) {
		ffts_[i].destroy();
		ffts_lt_[i].destroy();
	}
	delete[] ffts_;
	delete[] ffts_lt_;

	delete[] window_;
	delete[] coarse_aperiodicity_;
	delete[] coarse_frequency_axis_;
}

void D4C::prepareForD4c(int fs)
{
	fs_ = fs;
	fft_size_d4c_ =
		pow(2, 1 + static_cast<int>(log(4.0 * fs_ / world::kFloorF0D4C + 1) / world::kLog2));
	number_of_aperiodicities_ =
		static_cast<int>(MyMinDouble(world::kUpperLimit, fs_ / 2.0 -
									 world::kFrequencyInterval) / world::kFrequencyInterval);

	// window function
	window_length_ = static_cast<int>(world::kFrequencyInterval * fft_size_d4c_ / fs_) * 2 + 1;
    window_ =  new double[window_length_];
	NuttallWindow(window_length_, window_);
	
#ifdef _OPENMP
	num_thread_ = omp_get_num_procs();
#else
	num_thread_ = 1;
#endif

	coarse_aperiodicity_ = new double[(number_of_aperiodicities_ + 2) * num_thread_];
	int head_ind;
	for (int i = 0; i < num_thread_; ++i) {
		head_ind = (number_of_aperiodicities_ + 2) * i;
		coarse_aperiodicity_[head_ind] = -60.0;
		coarse_aperiodicity_[head_ind + number_of_aperiodicities_ + 1] = - world::kMySafeGuardMinimum;
	}
	
	coarse_frequency_axis_ = new double[number_of_aperiodicities_ + 2];
	coarse_frequency_axis_[number_of_aperiodicities_ + 1] = fs_ / 2.0;
	for (int i = 0; i <= number_of_aperiodicities_; ++i) {
		coarse_frequency_axis_[i] = i * world::kFrequencyInterval;
	}

	ffts_ = new ForwardRealFFT[num_thread_];
	for (int i = 0; i < num_thread_; ++i) {
		ForwardRealFFT forward_real_fft;
		forward_real_fft.initialize(fft_size_d4c_);
		ffts_[i] = forward_real_fft;
	}

	// FFTs for Love Train
	lowest_f0_ = 40.0;
	fft_size_lt_ = pow(2, 1 + static_cast<int>(log(3.0 * fs_ / lowest_f0_ + 1) / world::kLog2));

	ffts_lt_ = new ForwardRealFFT[num_thread_];
	for (int i = 0; i < num_thread_; ++i) {
		ForwardRealFFT forward_real_fft;
		forward_real_fft.initialize(fft_size_lt_);
		ffts_lt_[i] = forward_real_fft;
	}	
}

void D4C::compute(
	const double *x, int x_length,
	const double *temporal_positions, const double *f0, int f0_length,
	int fft_size, double **aperiodicity)
{
	x_ = x;
	x_length_ = x_length;
	time_pos_ = temporal_positions;
	f0_ = f0;
	f0_length_ = f0_length;
	fft_size_ = fft_size;
	int bin_size = fft_size_ / 2 + 1;

	// initialize
	constexpr double init_val = 1.0 - world::kMySafeGuardMinimum;
	for (int i = 0; i < f0_length_; ++i) {
		for (int j = 0; j < bin_size; ++j) {
			aperiodicity[i][j] = init_val;
		}
	}
  
	// D4C Love Train (Aperiodicity of 0 Hz is given by the different algorithm)
	double *aperiodicity0 = new double[f0_length_];
	loveTrain(aperiodicity0);

	double *frequency_axis = new double[bin_size];
	for (int i = 0; i < bin_size; ++i) {
		frequency_axis[i] = static_cast<double>(i) * fs_ / fft_size_;
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_thread_)
#endif
	for (int i = 0; i < f0_length_; ++i) {
		if (f0_[i] == 0 || aperiodicity0[i] <= option_.threshold) { continue; }

#ifdef _OPENMP
		int thread_id = omp_get_thread_num();
#else
		constexpr int thread_id = 0;
#endif

		double* ca_head = coarse_aperiodicity_ + (number_of_aperiodicities_ + 2) * thread_id;
		
		generalBody(MyMaxDouble(world::kFloorF0D4C, f0[i]), temporal_positions[i],
					ffts_[thread_id], ca_head + 1);
		
		// Linear interpolation to convert the coarse aperiodicity into its
		// spectral representation. (getAperiodicity)
		double *tmp_aperiodicity = aperiodicity[i];
		interp1(coarse_frequency_axis_, ca_head,
				number_of_aperiodicities_ + 2, frequency_axis, bin_size,
				tmp_aperiodicity);

		for_each(tmp_aperiodicity, tmp_aperiodicity + bin_size,
				 [](double& v){v = pow(10.0, v / 20.0);});
   	}
	
	delete[] aperiodicity0;
	delete[] frequency_axis;
}

//-----------------------------------------------------------------------------
// D4C::loveTrain() determines the aperiodicity with VUV detection.
// If a frame was determined as the unvoiced section, aperiodicity is set to
// very high value as the safeguard.
// If it was voiced section, the aperiodicity of 0 Hz is set to -60 dB.
//-----------------------------------------------------------------------------
void D4C::loveTrain(double *aperiodicity0)
{
	// Cumulative powers at 100, 4000, 7900 Hz are used for VUV identification.
	int boundary0 = static_cast<int>(ceil(100.0 * fft_size_lt_ / fs_));
	int boundary1 = static_cast<int>(ceil(4000.0 * fft_size_lt_ / fs_));
	int boundary2 = static_cast<int>(ceil(7900.0 * fft_size_lt_ / fs_));

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_thread_)
#endif
	for (int i = 0; i < f0_length_; ++i) {
		if (f0_[i] == 0.0) {
			aperiodicity0[i] = 0.0;
			continue;
		}

#ifdef _OPENMP
		int thread_id = omp_get_thread_num();
#else
		constexpr int thread_id = 0;
#endif
    
		aperiodicity0[i] = loveTrainSub(MyMaxDouble(f0_[i], lowest_f0_), time_pos_[i],
										fft_size_lt_, boundary0, boundary1, boundary2,
										ffts_lt_[thread_id]);
	}
}

double D4C::loveTrainSub(
	const double current_f0, const double current_position, const int fft_size,
	const int boundary0, const int boundary1, const int boundary2,
	const ForwardRealFFT &forward_real_fft)
{
	double *power_spectrum = new double[fft_size];

	int window_length = matlab_round(1.5 * fs_ / current_f0) * 2 + 1;
	getWindowedWaveform(current_f0, current_position,
						world::kBlackman, 3.0, forward_real_fft.waveform);

	for_each(forward_real_fft.waveform + window_length,
			 forward_real_fft.waveform + fft_size, [](double &v){v = 0;});
  
	fft_execute(forward_real_fft.forward_fft);
  
	for_each(power_spectrum, power_spectrum + boundary0 + 1, [](double &v){v = 0;});
  
	for (int i = boundary0 + 1; i < fft_size / 2 + 1; ++i) {
		power_spectrum[i] =
			forward_real_fft.spectrum[i][0] * forward_real_fft.spectrum[i][0] +
			forward_real_fft.spectrum[i][1] * forward_real_fft.spectrum[i][1];
	}
  
	for (int i = boundary0; i <= boundary2; ++i) { power_spectrum[i] += power_spectrum[i - 1]; }
  
	double aperiodicity0 = power_spectrum[boundary1] / power_spectrum[boundary2];
  
	delete[] power_spectrum;

	return aperiodicity0;
}

//-----------------------------------------------------------------------------
// D4C::getWindowedWaveform() windows the waveform by F0-adaptive window
// In the variable window_type, 1: hanning, 2: blackman
//-----------------------------------------------------------------------------
void D4C::getWindowedWaveform(
	double current_f0, double current_position, int window_type,
	double window_length_ratio, double *waveform)
{
	int half_window_length = matlab_round(window_length_ratio * fs_ / current_f0 / 2.0);

	int *base_index = new int[half_window_length * 2 + 1];
	int *safe_index = new int[half_window_length * 2 + 1];
	double *window  = new double[half_window_length * 2 + 1];

	// set parameters
	int n = - half_window_length;
	generate(base_index, base_index + half_window_length * 2 + 1, [&]{ return n++; });
  
	int origin = matlab_round(current_position * fs_ + 0.001);
	for (int i = 0; i <= half_window_length * 2; ++i)
    { safe_index[i] = MyMinInt(x_length_ - 1, MyMaxInt(0, origin + base_index[i])); }

	// Designing of the window function
	double position;
	double const_val1, const_val2;
	const_val1 = 2.0 / window_length_ratio / fs_;
	// Hanning window
	if (window_type == world::kHanning) {
		const_val2 = world::kPi * current_f0;
		for (int i = 0; i <= half_window_length * 2; ++i) {
			position = const_val1 * base_index[i];
			window[i] = 0.5 * cos(const_val2 * position) + 0.5;
		}
	}
	// Blackman window
	else {
		const_val2 = world::kPi * current_f0;
		for (int i = 0; i <= half_window_length * 2; ++i) {
			position = const_val1 * base_index[i];
			window[i] = 0.42 + 0.5 * cos(const_val2 * position) + 0.08 * cos(const_val2 * position * 2);
		}
	}
  
	// SetParametersForGetWindowedWaveform(half_window_length, x_length,
	//     current_position, fs, current_f0, window_type, window_length_ratio,
	//     base_index, safe_index, window);

	// F0-adaptive windowing
	for (int i = 0; i <= half_window_length * 2; ++i)
    { waveform[i] = x_[safe_index[i]] * window[i] + randn() * world::kMySafeGuardMinimum; }

	double tmp_weight1 = accumulate(waveform, waveform + half_window_length * 2 + 1, 0.0);
	double tmp_weight2 = accumulate(window, window + half_window_length * 2 + 1, 0.0);
	double weighting_coefficient = tmp_weight1 / tmp_weight2;
  
	for (int i = 0; i <= half_window_length * 2; ++i)
    { waveform[i] -= window[i] * weighting_coefficient; }

	delete[] base_index;
	delete[] safe_index;
	delete[] window;
}

//-----------------------------------------------------------------------------
// D4C::generalBody() calculates a spectral envelope at a temporal position.
//-----------------------------------------------------------------------------
void D4C::generalBody(
	const double current_f0, const double current_position,
	const ForwardRealFFT &forward_real_fft,
	double *coarse_aperiodicity)
{
	double *static_centroid = new double[fft_size_d4c_ / 2 + 1];
	double *smoothed_power_spectrum = new double[fft_size_d4c_ / 2 + 1];
	double *static_group_delay = new double[fft_size_d4c_ / 2 + 1];
  
	getStaticCentroid(current_f0, current_position, forward_real_fft, static_centroid);

	getSmoothedPowerSpectrum(current_f0, current_position, forward_real_fft, smoothed_power_spectrum);

   	getStaticGroupDelay(static_centroid, smoothed_power_spectrum, current_f0, static_group_delay);

   	getCoarseAperiodicity(static_group_delay, forward_real_fft, coarse_aperiodicity);
	
	// Revision of the result based on the F0
	double const_val = (current_f0 - 100) / 50.0;
	for (int i = 0; i < number_of_aperiodicities_; ++i)
    { coarse_aperiodicity[i] = MyMinDouble(0.0, coarse_aperiodicity[i] + const_val); }

	delete[] static_centroid;
	delete[] smoothed_power_spectrum;
	delete[] static_group_delay;
}

//-----------------------------------------------------------------------------
// GetStaticCentroid() calculates the temporally static energy centroid.
// Basic idea was proposed by H. Kawahara.
//-----------------------------------------------------------------------------
void D4C::getStaticCentroid(
	double current_f0, double current_position,
	const ForwardRealFFT &forward_real_fft, double *static_centroid)
{
	int bin_size = fft_size_d4c_ / 2 + 1;
  
	double *centroid1 = new double[bin_size];
	double *centroid2 = new double[bin_size];

	getCentroid(current_f0, current_position - 0.25 / current_f0,
				forward_real_fft, centroid1);
	getCentroid(current_f0, current_position + 0.25 / current_f0,
				forward_real_fft, centroid2);

	copy(centroid1, centroid1 + bin_size, static_centroid);
	for (int i = 0; i < bin_size; ++i) { static_centroid[i] += centroid2[i]; }

	DCCorrection(static_centroid, current_f0, fs_, fft_size_d4c_, static_centroid);
  
	delete[] centroid1;
	delete[] centroid2;
}

//-----------------------------------------------------------------------------
// D4C::getCentroid() calculates the energy centroid (see the book, time-frequency
// analysis written by L. Cohen).
//-----------------------------------------------------------------------------
void D4C::getCentroid(
	const double current_f0, const double current_position,
	const ForwardRealFFT &forward_real_fft, double *centroid)
{
	double *waveform = forward_real_fft.waveform;
	fft_complex *spectrum = forward_real_fft.spectrum;

	for_each(waveform, waveform + fft_size_d4c_, [](double &v){v = 0;});
  
	getWindowedWaveform(current_f0, current_position, world::kBlackman,
						4.0, waveform);
  
	double power = 0.0;
	for (int i = 0; i <= matlab_round(2.0 * fs_ / current_f0) * 2; ++i)
    { power += waveform[i] * waveform[i]; }

	power = sqrt(power);
	for (int i = 0; i <= matlab_round(2.0 * fs_ / current_f0) * 2; ++i)
    { waveform[i] /= power; }

	fft_execute(forward_real_fft.forward_fft);

	int bin_size = fft_size_d4c_ / 2 + 1;
	double *tmp_real = new double[bin_size];
	double *tmp_imag = new double[bin_size];
	for (int i = 0; i < bin_size; ++i) {
		tmp_real[i] = spectrum[i][0];
		tmp_imag[i] = spectrum[i][1];
	}

	for (int i = 0; i < fft_size_d4c_; ++i) { waveform[i] *= i + 1.0; }
  
	fft_execute(forward_real_fft.forward_fft);
  
	for (int i = 0; i < bin_size; ++i)
    { centroid[i] = spectrum[i][0] * tmp_real[i] + tmp_imag[i] * spectrum[i][1]; }

	delete[] tmp_real;
	delete[] tmp_imag;
}

//-----------------------------------------------------------------------------
// GetSmoothedPowerSpectrum() calculates the smoothed power spectrum.
// The parameters used for smoothing are optimized in davance.
//-----------------------------------------------------------------------------
void D4C::getSmoothedPowerSpectrum(
	const double current_f0, const double current_position,
	const ForwardRealFFT &forward_real_fft,
	double *smoothed_power_spectrum)
{
	double *waveform = forward_real_fft.waveform;
	fft_complex *spectrum = forward_real_fft.spectrum;

	for_each(waveform, waveform + fft_size_d4c_, [](double &v){v = 0;});
  
	getWindowedWaveform(current_f0, current_position, world::kHanning, 4.0, waveform);

	fft_execute(forward_real_fft.forward_fft);
  
	for (int i = 0; i <= fft_size_d4c_ / 2; ++i)
    { smoothed_power_spectrum[i] =
			spectrum[i][0] * spectrum[i][0] + spectrum[i][1] * spectrum[i][1]; }
  
	DCCorrection(smoothed_power_spectrum, current_f0, fs_, fft_size_d4c_,
				 smoothed_power_spectrum);
  
	LinearSmoothing(smoothed_power_spectrum, current_f0, fs_, fft_size_d4c_,
					smoothed_power_spectrum);
}

//-----------------------------------------------------------------------------
// D4C::getStaticGroupDelay() calculates the temporally static group delay.
// This is the fundamental parameter in D4C.
//-----------------------------------------------------------------------------
void D4C::getStaticGroupDelay(
	const double *static_centroid, const double *smoothed_power_spectrum,
	const double current_f0, double *static_group_delay)
{
	int bin_size = fft_size_d4c_ / 2 + 1;
  
	for (int i = 0; i < bin_size; ++i)
    { static_group_delay[i] = static_centroid[i] / smoothed_power_spectrum[i]; }
  
	LinearSmoothing(static_group_delay, current_f0 / 2.0, fs_, fft_size_d4c_,
					static_group_delay);

	double *smoothed_group_delay = new double[bin_size];
	LinearSmoothing(static_group_delay, current_f0, fs_, fft_size_d4c_,
					smoothed_group_delay);

	for (int i = 0; i < bin_size; ++i)
    { static_group_delay[i] -= smoothed_group_delay[i]; }

	delete[] smoothed_group_delay;
}

//-----------------------------------------------------------------------------
// GetCoarseAperiodicity() calculates the aperiodicity in multiples of 3 kHz.
// The upper limit is given based on the sampling frequency.
//-----------------------------------------------------------------------------
void D4C::getCoarseAperiodicity(
	const double *static_group_delay,
	const ForwardRealFFT &forward_real_fft,
	double *coarse_aperiodicity)
{
	double *waveform = forward_real_fft.waveform;
	fft_complex *spectrum = forward_real_fft.spectrum;

	for_each(waveform, waveform + fft_size_d4c_, [](double &v){v = 0;});
  
	int boundary = matlab_round(fft_size_d4c_ * 8.0 / window_length_);
	int half_window_length = window_length_ / 2;
	
	int bin_size = fft_size_d4c_ / 2 + 1;
	double *power_spectrum = new double[bin_size];
	int center;
  
	for (int i = 0; i < number_of_aperiodicities_; ++i) {
		center = static_cast<int>(world::kFrequencyInterval * (i + 1) * fft_size_d4c_ / fs_);
    
		for (int j = 0; j <= half_window_length * 2; ++j)
		{ waveform[j] = static_group_delay[center - half_window_length + j] * window_[j]; }
    
		fft_execute(forward_real_fft.forward_fft);
    
		for (int j = 0 ; j < bin_size; ++j)
		{ power_spectrum[j] = spectrum[j][0] * spectrum[j][0] + spectrum[j][1] * spectrum[j][1]; }
    
		std::sort(power_spectrum, power_spectrum + bin_size);
    
		for (int j = 1 ; j < bin_size; ++j) { power_spectrum[j] += power_spectrum[j - 1]; }
    
		coarse_aperiodicity[i] = 10 * log10(power_spectrum[bin_size - boundary - 2]
											/ power_spectrum[bin_size - 1]);
	}
  
	delete[] power_spectrum;
}


}
