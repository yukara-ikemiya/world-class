//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise,
// Copyright 2018 Yukara Ikemiya
//
// Voice synthesis based on f0, spectrogram and aperiodicity.
//-----------------------------------------------------------------------------

#include "synthesis.hpp"

#include <cmath>
#include <numeric>
#include <algorithm>
#include <omp.h>

#include "world_common.hpp"
#include "world_constantnumbers.hpp"
#include "world_matlabfunctions.hpp"


namespace
{
using std::accumulate;
using std::for_each;
using std::max_element;
}

namespace world_class
{

Synthesis::Synthesis(int fs, int fft_size, double frame_period)
	: fs_(fs), fft_size_(fft_size), frame_period_(frame_period / 1000.)
	, dc_remover_(new double[fft_size_])
{
#ifdef _OPENMP
	num_thread_ = omp_get_num_procs();
#else
	num_thread_ = 1;
#endif

	minimum_phase_ = new MinimumPhaseAnalysis[num_thread_];
	inverse_real_fft_ = new InverseRealFFT[num_thread_];
	forward_real_fft_ = new ForwardRealFFT[num_thread_];

	for (int ii = 0; ii < num_thread_; ii++) {
		minimum_phase_[ii].initialize(fft_size_);
		inverse_real_fft_[ii].initialize(fft_size_);
		forward_real_fft_[ii].initialize(fft_size_);
	}
	
	spectral_envelope_ = new double[(fft_size_ / 2 + 1) * num_thread_];
	aperiodic_ratio_ = new double[(fft_size_ / 2 + 1) * num_thread_];
	periodic_response_ = new double[fft_size_ * num_thread_];
	aperiodic_response_ = new double[fft_size_ * num_thread_];
	
	getDCRemover(fft_size_, dc_remover_);
}

Synthesis::~Synthesis()
{
	delete[] dc_remover_;
	delete[] spectral_envelope_;
	delete[] aperiodic_ratio_;
	delete[] periodic_response_;
	delete[] aperiodic_response_;

	for (int ii = 0; ii < num_thread_; ii++) {
		minimum_phase_[ii].destroy();
		inverse_real_fft_[ii].destroy();
		forward_real_fft_[ii].destroy();
	}

	delete[] minimum_phase_;
	delete[] inverse_real_fft_;
	delete[] forward_real_fft_;
}

void Synthesis::compute(
	const double *f0, int f0_length,
	const double * const *spectrogram, const double * const *aperiodicity,
	int out_length, double *out)
{
	// init
	for_each(out, out + out_length, [](double& v){v = 0;});

	double max_f0 = *max_element(f0, f0 + f0_length);
	
	int min_pulse_interval = (int)(fs_ / max_f0);
	int max_num_pulses = out_length / min_pulse_interval;

	double *pulse_locations = new double[max_num_pulses];
	int *pulse_locations_index = new int[max_num_pulses];
	double *pulse_locations_time_shift = new double[max_num_pulses];
	double *interpolated_vuv = new double[out_length];

	int number_of_pulses =
		getTimeBase(f0, f0_length, fs_, frame_period_,
					out_length, fs_ / fft_size_ + 1.0, pulse_locations, pulse_locations_index,
					pulse_locations_time_shift, interpolated_vuv);
	
#ifdef _OPENMP
	double *impulse_response = new double[fft_size_ * number_of_pulses];
#pragma omp parallel for num_threads(num_thread_)
	for (int ii = 0; ii < number_of_pulses; ii++) {
		int thread_id = omp_get_thread_num();
		
		int noise_size = pulse_locations_index[MyMinInt(number_of_pulses - 1, ii + 1)] -
						 pulse_locations_index[ii];
		
		getOneFrameSegment(
			interpolated_vuv[pulse_locations_index[ii]], noise_size,
			spectrogram, fft_size_, aperiodicity, f0_length, frame_period_,
			pulse_locations[ii], pulse_locations_time_shift[ii], fs_,
			forward_real_fft_ + thread_id, inverse_real_fft_ + thread_id, minimum_phase_ + thread_id,
			dc_remover_, impulse_response + ii * fft_size_, thread_id
		);
	}
	
	int index, b_index, e_index, len_valid;
	double *head = impulse_response;
	for (int ii = 0; ii < number_of_pulses; ii++) {
		index = pulse_locations_index[ii] - fft_size_ / 2;
		
		if (index + fft_size_ < 0 || index + 1 >= out_length) {
			head += fft_size_;
			continue;
		}
		
		b_index = (index + 1 < 0) ? abs(index + 1) : 0;
		e_index = (index + fft_size_ >= out_length) ? out_length - index - 1 : fft_size_;
		len_valid = e_index - b_index;
		index += b_index;
		
		for (int jj = b_index; jj < b_index + len_valid; jj++) {
			index++;
			out[index] += head[jj];
		}
		
		head += fft_size_;
	}	
#else
	double *impulse_response = new double[fft_size_];
	int noise_size;
	int index, b_index, e_index, len_valid;
	for (int ii = 0; ii < number_of_pulses; ii++) {
		noise_size = pulse_locations_index[MyMinInt(number_of_pulses - 1, ii + 1)] -
					 pulse_locations_index[ii];
		
		getOneFrameSegment(
			interpolated_vuv[pulse_locations_index[ii]], noise_size,
			spectrogram, fft_size_, aperiodicity, f0_length, frame_period_,
			pulse_locations[ii], pulse_locations_time_shift[ii], fs_,
			forward_real_fft_, inverse_real_fft_, minimum_phase_,
			dc_remover_, impulse_response
		);

		index = pulse_locations_index[ii] - fft_size_ / 2;

		if (index + fft_size_ < 0 || index + 1 >= out_length) { continue; }

		b_index = (index + 1 < 0) ? abs(index + 1) : 0;
		e_index = (index + fft_size_ >= out_length) ? out_length - index - 1 : fft_size_;
		len_valid = e_index - b_index;
		index += b_index;
		
		for (int jj = b_index; jj < b_index + len_valid; jj++) {
			index++;
			out[index] += impulse_response[jj];
		}
	}
#endif
	
	delete[] pulse_locations;
	delete[] pulse_locations_index;
	delete[] pulse_locations_time_shift;
	delete[] interpolated_vuv;
	delete[] impulse_response;
}


int Synthesis::getTimeBase(
	const double *f0, int f0_length, int fs,
	double frame_period, int y_length, double lowest_f0,
	double *pulse_locations, int *pulse_locations_index,
	double *pulse_locations_time_shift, double *interpolated_vuv)
{
	double *time_axis = new double[y_length];
	double *coarse_time_axis = new double[f0_length + 1];
	double *coarse_f0 = new double[f0_length + 1];
	double *coarse_vuv = new double[f0_length + 1];
	double *interpolated_f0 = new double[y_length];
	
	getTemporalParametersForTimeBase(
		f0, f0_length, fs, y_length, frame_period,
		lowest_f0, time_axis, coarse_time_axis, coarse_f0, coarse_vuv
	);
	
	interp1(coarse_time_axis, coarse_f0, f0_length + 1,
			time_axis, y_length, interpolated_f0);
	
	interp1(coarse_time_axis, coarse_vuv, f0_length + 1,
			time_axis, y_length, interpolated_vuv);
	
	for (int ii = 0; ii < y_length; ii++) {
		interpolated_vuv[ii] = interpolated_vuv[ii] > 0.5 ? 1.0 : 0.0;
		interpolated_f0[ii] =
			interpolated_vuv[ii] == 0.0 ? world::kDefaultF0 : interpolated_f0[ii];
	}
	
	int number_of_pulses =
		getPulseLocationsForTimeBase(
			interpolated_f0,
			time_axis, y_length, fs, pulse_locations, pulse_locations_index,
			pulse_locations_time_shift
		);
	
	delete[] coarse_vuv;
	delete[] coarse_f0;
	delete[] coarse_time_axis;
	delete[] time_axis;
	delete[] interpolated_f0;

	return number_of_pulses;
}

void Synthesis::getTemporalParametersForTimeBase(
	const double *f0, int f0_length,
    int fs, int y_length, double frame_period, double lowest_f0,
    double *time_axis, double *coarse_time_axis, double *coarse_f0,
    double *coarse_vuv)
{
	for (int ii = 0; ii < y_length; ii++) { time_axis[ii] = ii / (double)fs; }
	
	// the array 'coarse_time_axis' is supposed to have 'f0_length + 1' positions
	for (int ii = 0; ii < f0_length; ii++) {
		coarse_time_axis[ii] = ii * frame_period;
		coarse_f0[ii] = (f0[ii] < lowest_f0) ? 0.0 : f0[ii];
		coarse_vuv[ii] = (coarse_f0[ii] == 0.0) ? 0.0 : 1.0;
	}
	
	coarse_time_axis[f0_length] = f0_length * frame_period;
	coarse_f0[f0_length] = coarse_f0[f0_length - 1] * 2 - coarse_f0[f0_length - 2];
	coarse_vuv[f0_length] = coarse_vuv[f0_length - 1] * 2 - coarse_vuv[f0_length - 2];
}

int Synthesis::getPulseLocationsForTimeBase(
	const double *interpolated_f0,
    const double *time_axis, int y_length, int fs, double *pulse_locations,
    int *pulse_locations_index, double *pulse_locations_time_shift)
{
	double *total_phase = new double[y_length];
	double *wrap_phase = new double[y_length];
	double *wrap_phase_abs = new double[y_length - 1];

	double two_pi = 2.0 * world::kPi;
	double const_val = two_pi / fs;
	
	total_phase[0] = interpolated_f0[0] * const_val;
	wrap_phase[0] = fmod(total_phase[0], two_pi);
	
	for (int ii = 1; ii < y_length; ii++) {
		total_phase[ii] = total_phase[ii - 1] + interpolated_f0[ii] * const_val;
		wrap_phase[ii] = fmod(total_phase[ii], two_pi);
		wrap_phase_abs[ii - 1] = fabs(wrap_phase[ii] - wrap_phase[ii - 1]);
	}

	int number_of_pulses = 0;
	double y1, y2, x;
	for (int ii = 0; ii < y_length - 1; ii++) {
		if (wrap_phase_abs[ii] > world::kPi) {
			pulse_locations[number_of_pulses] = time_axis[ii];
			pulse_locations_index[number_of_pulses] = ii;

			// calculate an exact phase of a zero cross point
			y1 = wrap_phase[ii] - two_pi;
			y2 = wrap_phase[ii + 1];
			x = - y1 / (y2 - y1);
			pulse_locations_time_shift[number_of_pulses] = x / fs;

			++number_of_pulses;
		}
	}

	delete[] wrap_phase_abs;
	delete[] wrap_phase;
	delete[] total_phase;

	return number_of_pulses;
}

void Synthesis::getDCRemover(int fft_size, double *dc_remover)
{
	double const_val = 2.0 * world::kPi / (1.0 + fft_size);
	for (int ii = 0; ii < fft_size / 2; ii++) {
		dc_remover[ii] = 0.5 - 0.5 * cos(const_val * (ii + 1.0));
	}

	double dc_component = accumulate(dc_remover, dc_remover + fft_size / 2, 0.0) * 2;
	
	for (int ii = 0; ii < fft_size / 2; ii++) {
		dc_remover[ii] /= dc_component;
		dc_remover[fft_size - ii - 1] = dc_remover[ii];
	}
}

//-----------------------------------------------------------------------------
// GetOneFrameSegment() calculates a periodic and aperiodic response at a time.
//-----------------------------------------------------------------------------
void Synthesis::getOneFrameSegment(
	double current_vuv, int noise_size,
	const double * const *spectrogram, int fft_size,
	const double * const *aperiodicity, int f0_length, double frame_period,
	double current_time, double fractional_time_shift, int fs,
	const ForwardRealFFT *forward_real_fft, const InverseRealFFT *inverse_real_fft,
	MinimumPhaseAnalysis *minimum_phase, const double *dc_remover,
	double *response, int thread_id)
{
	double *spectral_envelope = spectral_envelope_ + (fft_size_ / 2 + 1) * thread_id;
	double *aperiodic_ratio = aperiodic_ratio_ + (fft_size_ / 2 + 1) * thread_id;

	double *periodic_response = periodic_response_ + fft_size_ * thread_id;
	double *aperiodic_response = aperiodic_response_ + fft_size_ * thread_id;
	
	getSpectralEnvelope(current_time, frame_period, f0_length, spectrogram,
						fft_size, spectral_envelope);
	
	getAperiodicRatio(current_time, frame_period, f0_length, aperiodicity,
					  fft_size, aperiodic_ratio);

	// synthesis of the periodic response
	getPeriodicResponse(fft_size, spectral_envelope, aperiodic_ratio,
						current_vuv, inverse_real_fft, minimum_phase, dc_remover,
						fractional_time_shift, fs, periodic_response);

	// synthesis of the aperiodic response
	getAperiodicResponse(noise_size, fft_size, spectral_envelope,
						 aperiodic_ratio, current_vuv, forward_real_fft,
						 inverse_real_fft, minimum_phase, aperiodic_response);

	double sqrt_noise_size = sqrt((double)noise_size);
	for (int ii = 0; ii < fft_size; ii++) {
		response[ii] =
			(periodic_response[ii] * sqrt_noise_size + aperiodic_response[ii]) / fft_size;
	}
}

void Synthesis::getSpectralEnvelope(
	double current_time, double frame_period,
    int f0_length, const double * const *spectrogram, int fft_size,
    double *spectral_envelope)
{
	int current_frame_floor = MyMinInt(f0_length - 1, (int)floor(current_time / frame_period));
	int current_frame_ceil = MyMinInt(f0_length - 1, (int)ceil(current_time / frame_period));
	double interpolation = current_time / frame_period - current_frame_floor;

	if (current_frame_floor == current_frame_ceil) {
		const double *floor_head = spectrogram[current_frame_floor];
		for (int ii = 0; ii <= fft_size / 2; ii++) {
			spectral_envelope[ii] = fabs(floor_head[ii]);
		}
	} else {
		const double *floor_head = spectrogram[current_frame_floor];
		const double *ceil_head = spectrogram[current_frame_ceil];
		for (int ii = 0; ii <= fft_size / 2; ii++) {
			spectral_envelope[ii] = 
				(1.0 - interpolation) * fabs(floor_head[ii]) + interpolation * fabs(ceil_head[ii]);
		}
	}
}

void Synthesis::getAperiodicRatio(
	double current_time, double frame_period,
    int f0_length, const double * const *aperiodicity, int fft_size,
    double *aperiodic_spectrum)
{
	int current_frame_floor = MyMinInt(f0_length - 1, (int)floor(current_time / frame_period));
	int current_frame_ceil = MyMinInt(f0_length - 1, (int)ceil(current_time / frame_period));
	double interpolation = current_time / frame_period - current_frame_floor;

	if (current_frame_floor == current_frame_ceil) {
		const double *floor_head = aperiodicity[current_frame_floor];
		for (int i = 0; i <= fft_size / 2; ++i) {
			aperiodic_spectrum[i] = pow(GetSafeAperiodicity(floor_head[i]), 2.0);
		}
	} else {
		const double *floor_head = aperiodicity[current_frame_floor];
		const double *ceil_head = aperiodicity[current_frame_ceil];
		for (int i = 0; i <= fft_size / 2; ++i) {
			aperiodic_spectrum[i] = pow(
				(1.0 - interpolation) * getSafeAperiodicity(floor_head[i]) +
				interpolation * getSafeAperiodicity(ceil_head[i]), 2.0);
		}
	}
}

inline double Synthesis::getSafeAperiodicity(double x)
{
	return MyMaxDouble(0.001, MyMinDouble(0.999999999999, x));
}

//-----------------------------------------------------------------------------
// getPeriodicResponse() calculates a periodic response.
//-----------------------------------------------------------------------------
void Synthesis::getPeriodicResponse(
	int fft_size, const double *spectrum,
    const double *aperiodic_ratio, double current_vuv,
    const InverseRealFFT *inverse_real_fft,
    MinimumPhaseAnalysis *minimum_phase, const double *dc_remover,
    double fractional_time_shift, int fs, double *periodic_response)
{
	if (current_vuv <= 0.5 || aperiodic_ratio[0] > 0.999) {
		for_each(periodic_response, periodic_response + fft_size, [](double &v){v = 0;});
		return;
	}

	for (int ii = 0; ii <= fft_size / 2; ii++) {
		minimum_phase->log_spectrum[ii] = log(spectrum[ii] * (1.0 - aperiodic_ratio[ii]) +
											  world::kMySafeGuardMinimum) / 2.0;
	}
	
	minimum_phase->compute();

	for (int ii = 0; ii <= fft_size / 2; ii++) {
		inverse_real_fft->spectrum[ii][0] = minimum_phase->minimum_phase_spectrum[ii][0];
		inverse_real_fft->spectrum[ii][1] = minimum_phase->minimum_phase_spectrum[ii][1];
	}

	// apply fractional time delay of fractional_time_shift seconds
	// using linear phase shift
	double coefficient = 2.0 * world::kPi * fractional_time_shift * fs / fft_size;
	getSpectrumWithFractionalTimeShift(fft_size, coefficient, inverse_real_fft);

	fft_execute(inverse_real_fft->inverse_fft);
	
	fftshift(inverse_real_fft->waveform, fft_size, periodic_response);
	
	removeDCComponent(periodic_response, fft_size, dc_remover, periodic_response);
}

//-----------------------------------------------------------------------------
// getSpectrumWithFractionalTimeShift() calculates a periodic spectrum with
// the fractional time shift under 1/fs.
//-----------------------------------------------------------------------------
void Synthesis::getSpectrumWithFractionalTimeShift(
	int fft_size, double coefficient, const InverseRealFFT *inverse_real_fft)
{
	double re, im, re2, im2;
	for (int ii = 0; ii <= fft_size / 2; ii++) {
		fft_complex &tmp_sp = inverse_real_fft->spectrum[ii];
		re = tmp_sp[0];
		im = tmp_sp[1];
		re2 = cos(coefficient * ii);
		im2 = sqrt(1.0 - re2 * re2);  // sin(pshift)

		tmp_sp[0] = re * re2 - im * im2;
		tmp_sp[1] = re * im2 + im * re2;
	}
}

void Synthesis::removeDCComponent(
	const double *periodic_response, int fft_size,
    const double *dc_remover, double *new_periodic_response)
{
	int half_fft_size = fft_size / 2;
	
	double dc_component =
		accumulate(periodic_response + half_fft_size, periodic_response + fft_size, 0.0);

	double dc_remove_val;
	for (int ii = 0; ii < half_fft_size; ii++) {
		dc_remove_val = - dc_component * dc_remover[ii];
		new_periodic_response[ii] = dc_remove_val;
		new_periodic_response[ii + half_fft_size] += dc_remove_val;
	}
}

//-----------------------------------------------------------------------------
// getAperiodicResponse() calculates an aperiodic response.
//-----------------------------------------------------------------------------
void Synthesis::getAperiodicResponse(
	int noise_size, int fft_size,
    const double *spectrum, const double *aperiodic_ratio, double current_vuv,
    const ForwardRealFFT *forward_real_fft,
    const InverseRealFFT *inverse_real_fft,
    MinimumPhaseAnalysis *minimum_phase, double *aperiodic_response)
{
	getNoiseSpectrum(noise_size, fft_size, forward_real_fft);

	if (current_vuv != 0.0) {
		for (int ii = 0; ii <= fft_size / 2; ii++) {
			minimum_phase->log_spectrum[ii] = log(spectrum[ii] * aperiodic_ratio[ii]) / 2.0;
		}
	} else {
		for (int ii = 0; ii <= fft_size / 2; ii++) {
			minimum_phase->log_spectrum[ii] = log(spectrum[ii]) / 2.0;
		}
	}
	
	minimum_phase->compute();

	for (int ii = 0; ii <= fft_size / 2; ii++) {
		fft_complex &tmp_noise_spec = minimum_phase->minimum_phase_spectrum[ii];
		fft_complex &tmp_sp = forward_real_fft->spectrum[ii];
		fft_complex &tmp_ap = inverse_real_fft->spectrum[ii];
		
		tmp_ap[0] = tmp_noise_spec[0] * tmp_sp[0] - tmp_noise_spec[1] * tmp_sp[1];	
		tmp_ap[1] = tmp_noise_spec[0] * tmp_sp[1] + tmp_noise_spec[1] * tmp_sp[0];
	}
	
	fft_execute(inverse_real_fft->inverse_fft);
	
	fftshift(inverse_real_fft->waveform, fft_size, aperiodic_response);
}

void Synthesis::getNoiseSpectrum(
	int noise_size, int fft_size, const ForwardRealFFT *forward_real_fft)
{
	double *waveform = forward_real_fft->waveform;
	
	for (int ii = 0; ii < noise_size; ii++) {
		waveform[ii] = randn();
	}

	double average = accumulate(waveform, waveform + noise_size, 0.0);
	average /= noise_size;

	for_each(waveform, waveform + noise_size, [&](double &v){v -= average;});
	for_each(waveform + noise_size, waveform + fft_size, [](double &v){v = 0.0;});
	
	fft_execute(forward_real_fft->forward_fft);
}

}
