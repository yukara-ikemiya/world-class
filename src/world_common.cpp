//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise,
// Copyright 2018 Yukara Ikemiya
//
// world_common.cpp includes functions used in at least two files.
// (1) Common functions
// (2) FFT, IFFT and minimum phase analysis.
//
// In FFT analysis and minimum phase analysis,
// Functions "Initialize*()" allocate the mamory.
// Functions "Destroy*()" free the accolated memory.
// FFT size is used for initialization, and structs are used to keep the memory.
// Functions "GetMinimumPhaseSpectrum()" calculate minimum phase spectrum.
// Forward and inverse FFT do not have the function "Get*()",
// because forward FFT and inverse FFT can run in one step.
//
//-----------------------------------------------------------------------------
#include "world_common.hpp"

#include <cmath>

#include "world_constantnumbers.hpp"
#include "world_matlabfunctions.hpp"

namespace
{
static void SetParametersForLinearSmoothing(
	int boundary, int fft_size, int fs,
    double width, const double *power_spectrum, double *mirroring_spectrum,
    double *mirroring_segment, double *frequency_axis)
{
	for (int i = 0; i < boundary; ++i)
	{ mirroring_spectrum[i] = power_spectrum[boundary - i]; }
	
	for (int i = boundary; i < fft_size / 2 + boundary; ++i)
	{ mirroring_spectrum[i] = power_spectrum[i - boundary]; }

	int st = fft_size / 2 + boundary;
	for (int i = st; i <= st + boundary; ++i) {
		mirroring_spectrum[i] =
			power_spectrum[fft_size / 2 - (i - (fft_size / 2 + boundary))];
	}

	mirroring_segment[0] = mirroring_spectrum[0] * fs / fft_size;
	for (int i = 1; i < fft_size / 2 + boundary * 2 + 1; ++i) {
		mirroring_segment[i] = mirroring_spectrum[i] * fs / fft_size +
							   mirroring_segment[i - 1];
	}

	for (int i = 0; i <= fft_size / 2; ++i)
	{ frequency_axis[i] = static_cast<double>(i) / fft_size * fs - width / 2.0; }
}
}  // namespace

//-----------------------------------------------------------------------------
int GetSuitableFFTSize(int sample) {
	return static_cast<int>(pow(2.0,
								static_cast<int>(log(static_cast<double>(sample)) / world::kLog2) + 1.0));
}

void DCCorrection(const double *input, double f0, int fs, int fft_size,
				  double *output) {
	int upper_limit = 2 + static_cast<int>(f0 * fft_size / fs);
	double *low_frequency_replica = new double[upper_limit];
	double *low_frequency_axis = new double[upper_limit];

	for (int i = 0; i < upper_limit; ++i)
		low_frequency_axis[i] = static_cast<double>(i) * fs / fft_size;

	int upper_limit_replica = upper_limit - 1;
	interp1Q(f0 - low_frequency_axis[0],
			 -static_cast<double>(fs) / fft_size, input, upper_limit + 1,
			 low_frequency_axis, upper_limit_replica, low_frequency_replica);

	for (int i = 0; i < upper_limit_replica; ++i)
		output[i] = input[i] + low_frequency_replica[i];

	delete[] low_frequency_replica;
	delete[] low_frequency_axis;
}

void LinearSmoothing(const double *input, double width, int fs, int fft_size,
					 double *output) {
	int boundary = static_cast<int>(width * fft_size / fs) + 1;

	// These parameters are set by the other function.
	double *mirroring_spectrum = new double[fft_size / 2 + boundary * 2 + 1];
	double *mirroring_segment = new double[fft_size / 2 + boundary * 2 + 1];
	double *frequency_axis = new double[fft_size / 2 + 1];
	SetParametersForLinearSmoothing(boundary, fft_size, fs, width,
									input, mirroring_spectrum, mirroring_segment, frequency_axis);

	double *low_levels = new double[fft_size / 2 + 1];
	double *high_levels = new double[fft_size / 2 + 1];
	double origin_of_mirroring_axis = -(boundary - 0.5) * fs / fft_size;
	double discrete_frequency_interval = static_cast<double>(fs) / fft_size;

	interp1Q(origin_of_mirroring_axis, discrete_frequency_interval,
			 mirroring_segment, fft_size / 2 + boundary * 2 + 1, frequency_axis,
			 fft_size / 2 + 1, low_levels);

	for (int i = 0; i <= fft_size / 2; ++i) frequency_axis[i] += width;

	interp1Q(origin_of_mirroring_axis, discrete_frequency_interval,
			 mirroring_segment, fft_size / 2 + boundary * 2 + 1, frequency_axis,
			 fft_size / 2 + 1, high_levels);

	for (int i = 0; i <= fft_size / 2; ++i)
		output[i] = (high_levels[i] - low_levels[i]) / width;

	delete[] mirroring_spectrum;
	delete[] mirroring_segment;
	delete[] frequency_axis;
	delete[] low_levels;
	delete[] high_levels;
}

void NuttallWindow(int y_length, double *y) {
	double tmp;
	for (int i = 0; i < y_length; ++i) {
		tmp  = i / (y_length - 1.0);
		y[i] = 0.355768 - 0.487396 * cos(2.0 * world::kPi * tmp) +
			   0.144232 * cos(4.0 * world::kPi * tmp) -
			   0.012604 * cos(6.0 * world::kPi * tmp);
	}
}

//-----------------------------------------------------------------------------
// FFT, IFFT and minimum phase analysis

void ForwardRealFFT::initialize(const int n)
{
	fft_size = n;
	waveform = new double[fft_size];
	spectrum = new fft_complex[fft_size / 2 + 1];
	forward_fft = fft_plan_dft_r2c_1d(fft_size, waveform, spectrum, FFT_ESTIMATE);
}

void ForwardRealFFT::destroy()
{
	fft_destroy_plan(forward_fft);
	delete[] spectrum;
	delete[] waveform;
}

void InverseRealFFT::initialize(const int n)
{
	fft_size = n;
	waveform = new double[fft_size];
	spectrum = new fft_complex[fft_size / 2 + 1];
	inverse_fft = fft_plan_dft_c2r_1d(fft_size, spectrum, waveform, FFT_ESTIMATE);
}

void InverseRealFFT::destroy()
{
	fft_destroy_plan(inverse_fft);
	delete[] spectrum;
	delete[] waveform;
}

void InverseComplexFFT::initialize(const int n)
{
	fft_size = n;
	input = new fft_complex[fft_size];
	output = new fft_complex[fft_size];
	inverse_fft = fft_plan_dft_1d(fft_size, input, output, FFT_BACKWARD, FFT_ESTIMATE);
}

void InverseComplexFFT::destroy()
{
	fft_destroy_plan(inverse_fft);
	delete[] input;
	delete[] output;
}

void MinimumPhaseAnalysis::initialize(const int n)
{
	fft_size = n;
	log_spectrum = new double[fft_size];
	minimum_phase_spectrum = new fft_complex[fft_size];
	cepstrum = new fft_complex[fft_size];
	inverse_fft = fft_plan_dft_r2c_1d(fft_size, log_spectrum, cepstrum, FFT_ESTIMATE);
	forward_fft = fft_plan_dft_1d(fft_size, cepstrum, minimum_phase_spectrum,
								  FFT_FORWARD, FFT_ESTIMATE);
}

void MinimumPhaseAnalysis::destroy()
{
	fft_destroy_plan(forward_fft);
	fft_destroy_plan(inverse_fft);
	delete[] cepstrum;
	delete[] log_spectrum;
	delete[] minimum_phase_spectrum;
}

void MinimumPhaseAnalysis::compute()
{
	// Mirroring
	for (int i = fft_size / 2 + 1; i < fft_size; ++i)
	{ log_spectrum[i] = log_spectrum[fft_size - i]; }

	// This fft_plan carries out "forward" FFT.
	// To carriy out the Inverse FFT, the sign of imaginary part
	// is inverted after FFT.
	fft_execute(inverse_fft);
  
	cepstrum[0][1] *= -1.0;
	for (int i = 1; i < fft_size / 2; ++i) {
		cepstrum[i][0] *= 2.0;
		cepstrum[i][1] *= -2.0;
	}
  
	cepstrum[fft_size / 2][1] *= -1.0;
	for (int i = fft_size / 2 + 1; i < fft_size; ++i) {
		cepstrum[i][0] = 0.0;
		cepstrum[i][1] = 0.0;
	}

	fft_execute(forward_fft);

	// Since x is complex number, calculation of exp(x) is as following.
	// Note: This FFT library does not keep the aliasing.
	double tmp;
	for (int i = 0; i <= fft_size / 2; ++i) {
		tmp = exp(minimum_phase_spectrum[i][0] / fft_size);
		minimum_phase_spectrum[i][0] = tmp *
									   cos(minimum_phase_spectrum[i][1] /
										   fft_size);
		minimum_phase_spectrum[i][1] = tmp *
									   sin(minimum_phase_spectrum[i][1] /
										   fft_size);
	}
}

