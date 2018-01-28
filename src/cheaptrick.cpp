//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise,
// Copyright 2018 Yukara Ikemiya
//
// Spectral envelope estimation on the basis of the idea of CheapTrick.
//-----------------------------------------------------------------------------

#include "cheaptrick.hpp"

#include <cmath>
#include <algorithm>

#include "world_common.hpp"
#include "world_constantnumbers.hpp"
#include "world_matlabfunctions.hpp"

using namespace std;

namespace world_class
{

CheapTrickOption::CheapTrickOption()
	: q1(-0.15), f0_floor(world::kFloorF0), fft_size(0)
{}


CheapTrick::CheapTrick(const int fs)
{
	option_.fft_size = getFFTSizeForCheapTrick(fs, option_.f0_floor);

	fs_ = fs;
	f0_floor_ = getF0FloorForCheapTrick(fs_, option_.fft_size);
}

CheapTrick::CheapTrick(const int fs, const CheapTrickOption &option)
{
	option_.q1 = option.q1;
	option_.f0_floor = option.f0_floor;
	option_.fft_size = (option.fft_size == 0)
					   ? getFFTSizeForCheapTrick(fs, option_.f0_floor)
					   : option.fft_size;

	fs_ = fs;
	f0_floor_ = getF0FloorForCheapTrick(fs_, option_.fft_size);
}


void CheapTrick::compute(
	const double *x, int x_length, const double *temporal_positions,
	const double *f0, int f0_length, double **spectrogram)
{
	x_ = x;
	x_length_ = x_length;
	int bin_size = option_.fft_size / 2 + 1;
  
#ifndef _OPENMP
	double *spectral_envelope = new double[option_.fft_size];
	ForwardRealFFT forward_real_fft;
	forward_real_fft.initialize(option_.fft_size);
	InverseRealFFT inverse_real_fft;
	inverse_real_fft.initialize(option_.fft_size);
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i < f0_length; ++i) {
#ifdef _OPENMP
		double *spectral_envelope = new double[option_.fft_size];
		ForwardRealFFT forward_real_fft;
		forward_real_fft.initialize(option_.fft_size);
		InverseRealFFT inverse_real_fft;
		inverse_real_fft.initialize(option_.fft_size);
#endif
    
		double current_f0 = (f0[i] <= f0_floor_) ? world::kDefaultF0 : f0[i];
    
		generalBody(current_f0, temporal_positions[i], forward_real_fft, inverse_real_fft,
					spectral_envelope);

		copy(spectral_envelope, spectral_envelope + bin_size, spectrogram[i]);

#ifdef _OPENMP
		forward_real_fft.destroy();
		inverse_real_fft.destroy();
		delete[] spectral_envelope;
#endif
	}

#ifndef _OPENMP
	forward_real_fft.destroy();
	inverse_real_fft.destroy();
	delete[] spectral_envelope;
#endif
}

int CheapTrick::getFFTSizeForCheapTrick(const int fs, const double f0_floor)
{
	return pow(2, 1 + static_cast<int>(log(3.0 * fs / f0_floor + 1) / world::kLog2));
}

double CheapTrick::getF0FloorForCheapTrick(const int fs, const int fft_size)
{
	return 3 * fs / (fft_size - 3.0);
}


void CheapTrick::generalBody(
	const double current_f0, const double current_position,
	const ForwardRealFFT &forward_real_fft,
	const InverseRealFFT &inverse_real_fft,
	double *spectral_envelope)
{
	// F0-adaptive windowing
	getWindowedWaveform(current_f0, current_position, forward_real_fft);
  
	// Calculate power spectrum with DC correction
	// Note: The calculated power spectrum is stored in an array for waveform.
	// In this imprementation, power spectrum is transformed by FFT (NOT IFFT).
	// However, the same result is obtained.
	// This is tricky but important for simple implementation.
	getPowerSpectrum(current_f0, forward_real_fft);

	// Smoothing of the power (linear axis)
	// forward_real_fft.waveform is the power spectrum.
	LinearSmoothing(forward_real_fft.waveform, current_f0 * 2.0 / 3.0,
					fs_, option_.fft_size, forward_real_fft.waveform);
  
	// Add infinitesimal noise
	// This is a safeguard to avoid including zero in the spectrum.
	addInfinitesimalNoise(forward_real_fft.waveform, forward_real_fft.waveform);
  
	// Smoothing (log axis) and spectral recovery on the cepstrum domain.
	smoothingWithRecovery(current_f0, forward_real_fft, inverse_real_fft, spectral_envelope);
}

void CheapTrick::getWindowedWaveform(
	const double current_f0, const double currnet_position,
	const ForwardRealFFT &forward_real_fft)
{
	int half_window_length = matlab_round(1.5 * fs_ / current_f0);

	int *base_index = new int[half_window_length * 2 + 1];
	int *safe_index = new int[half_window_length * 2 + 1];
	double *window  = new double[half_window_length * 2 + 1];

	setParametersForGetWindowedWaveform(half_window_length, currnet_position, current_f0,
										base_index, safe_index, window);

	// F0-adaptive windowing
	double *waveform = forward_real_fft.waveform;
	for (int i = 0; i <= half_window_length * 2; ++i)
		waveform[i] = x_[safe_index[i]] * window[i] + randn() * 0.000000000000001;
	double tmp_weight1 = 0;
	double tmp_weight2 = 0;
	for (int i = 0; i <= half_window_length * 2; ++i) {
		tmp_weight1 += waveform[i];
		tmp_weight2 += window[i];
	}
	double weighting_coefficient = tmp_weight1 / tmp_weight2;
	for (int i = 0; i <= half_window_length * 2; ++i)
		waveform[i] -= window[i] * weighting_coefficient;

	delete[] base_index;
	delete[] safe_index;
	delete[] window;
}

void CheapTrick::setParametersForGetWindowedWaveform(
	const int half_window_length, const double currnet_position, const double current_f0,
	int *base_index, int *safe_index, double *window)
{
	int n = - half_window_length;
	generate(base_index, base_index + half_window_length * 2 + 1, [&]{return n++;});
  
	// for (int i = -half_window_length; i <= half_window_length; ++i)
	//   base_index[i + half_window_length] = i;
  
	int origin = matlab_round(currnet_position * fs_ + 0.001);
	for (int i = 0; i <= half_window_length * 2; ++i)
    { safe_index[i] = MyMinInt(x_length_ - 1, MyMaxInt(0, origin + base_index[i])); }

	// Designing of the window function
	double average = 0.0;
	double position;
	for (int i = 0; i <= half_window_length * 2; ++i) {
		position = base_index[i] / 1.5 / fs_;
		window[i] = 0.5 * cos(world::kPi * position * current_f0) + 0.5;
		average += window[i] * window[i];
	}
  
	average = sqrt(average);

	for_each(window, window + half_window_length * 2 + 1, [&](double &v){v /= average;});
	// for (int i = 0; i <= half_window_length * 2; ++i) window[i] /= average;
}

void CheapTrick::getPowerSpectrum(
	const double f0,
	const ForwardRealFFT &forward_real_fft)
{
	int half_window_length = matlab_round(1.5 * fs_ / f0);

	// FFT
	for (int i = half_window_length * 2 + 1; i < option_.fft_size; ++i)
		forward_real_fft.waveform[i] = 0.0;
	fft_execute(forward_real_fft.forward_fft);

	// Calculation of the power spectrum.
	double *power_spectrum = forward_real_fft.waveform;
	for (int i = 0; i <= option_.fft_size / 2; ++i)
		power_spectrum[i] =
			forward_real_fft.spectrum[i][0] * forward_real_fft.spectrum[i][0] +
			forward_real_fft.spectrum[i][1] * forward_real_fft.spectrum[i][1];
  
	// DC correction
	DCCorrection(power_spectrum, f0, fs_, option_.fft_size, power_spectrum);
}

void CheapTrick::addInfinitesimalNoise(
	const double *input_spectrum,
	double *output_spectrum)
{
	int bin_size = option_.fft_size / 2 + 1;
	copy(input_spectrum, input_spectrum + bin_size, output_spectrum);
	for (int i = 0; i < bin_size; ++i)
    { output_spectrum[i] += fabs(randn()) * world::kEps; }
}

void CheapTrick::smoothingWithRecovery(
	const double f0,
	const ForwardRealFFT &forward_real_fft,
	const InverseRealFFT &inverse_real_fft,
	double *spectral_envelope)
{
	int fft_size = option_.fft_size;
	double q1 = option_.q1;
  
	double *smoothing_lifter = new double[fft_size];
	double *compensation_lifter = new double[fft_size];

	smoothing_lifter[0] = 1.0;
	compensation_lifter[0] = (1.0 - 2.0 * q1) + 2.0 * q1;
	double quefrency;
  
	for (int i = 1; i <= forward_real_fft.fft_size / 2; ++i) {
		quefrency = static_cast<double>(i) / fs_;
		smoothing_lifter[i] = sin(world::kPi * f0 * quefrency) /
							  (world::kPi * f0 * quefrency);
		compensation_lifter[i] = (1.0 - 2.0 * q1) + 2.0 * q1 *
								 cos(2.0 * world::kPi * quefrency * f0);
	}

	double *waveform = forward_real_fft.waveform;
	for (int i = 0; i <= fft_size / 2; ++i)
    { waveform[i] = log(waveform[i]); }

	reverse_copy(waveform + 1, waveform + fft_size / 2, waveform + fft_size / 2 + 1);
  
	fft_execute(forward_real_fft.forward_fft);

	for (int i = 0; i <= fft_size / 2; ++i) {
		inverse_real_fft.spectrum[i][0] = forward_real_fft.spectrum[i][0] *
										  smoothing_lifter[i] *
										  compensation_lifter[i] / fft_size;
		inverse_real_fft.spectrum[i][1] = 0.0;
	}
  
	fft_execute(inverse_real_fft.inverse_fft);

	for (int i = 0; i <= fft_size / 2; ++i)
    { spectral_envelope[i] = exp(inverse_real_fft.waveform[i]); }
  
	delete[] smoothing_lifter;
	delete[] compensation_lifter;
}

} // end namespace world_class
