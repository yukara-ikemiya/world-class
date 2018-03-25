//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise,
// Copyright 2018 Yukara Ikemiya
//
// Voice synthesis based on f0, spectrogram and aperiodicity.
//-----------------------------------------------------------------------------

#ifndef WORLD_CLASS_SYNTHESIS_HPP
#define WORLD_CLASS_SYNTHESIS_HPP

#include "world_common.hpp"

namespace world_class
{

class Synthesis
{

public:

	//-----------------------------------------------------------------------------
	// Constructor
	//
	// Input:
	//   fs                   : Sampling frequency
	//   fft_size             : FFT size
	//   frame_period         : Temporal period used for the analysis
	//-----------------------------------------------------------------------------
	Synthesis(int fs, int fft_size, double frame_period);
	~Synthesis();

	//-----------------------------------------------------------------------------
	// compute() synthesize the voice based on f0, spectrogram and
	// aperiodicity (not excitation signal).
	//
	// Input:
	//   f0                   : f0 contour
	//   f0_length            : Length of f0
	//   spectrogram          : Spectrogram estimated by CheapTrick
	//   aperiodicity         : Aperiodicity spectrogram based on D4C
	//   out_length           : Length of the output signal (Memory of y has been
	//                          allocated in advance)
	// Output:
	//   out                  : Calculated speech
	//-----------------------------------------------------------------------------
	
	void compute(
		const double *f0, int f0_length, 
		const double * const *spectrogram, const double * const *aperiodicity, 
		int out_length, double *out
	);

private:

	// common parameters
	int fs_;
	int fft_size_;
	double frame_period_;
	double *dc_remover_;
	int num_thread_;
	
	// buffers for computation
	MinimumPhaseAnalysis *minimum_phase_;
	InverseRealFFT *inverse_real_fft_;
	ForwardRealFFT *forward_real_fft_;
	double *spectral_envelope_;
	double *aperiodic_ratio_;
	double *periodic_response_;
	double *aperiodic_response_;
	
	int getTimeBase(
		const double *f0, int f0_length, int fs,
		double frame_period, int y_length, double lowest_f0,
		double *pulse_locations, int *pulse_locations_index,
		double *pulse_locations_time_shift, double *interpolated_vuv
	);

	void getTemporalParametersForTimeBase(
		const double *f0, int f0_length,
		int fs, int y_length, double frame_period, double lowest_f0,
		double *time_axis, double *coarse_time_axis, double *coarse_f0,
		double *coarse_vuv
	);

	int getPulseLocationsForTimeBase(
		const double *interpolated_f0,
		const double *time_axis, int y_length, int fs, double *pulse_locations,
		int *pulse_locations_index, double *pulse_locations_time_shift
	);

	void getDCRemover(int fft_size, double *dc_remover);

	void getOneFrameSegment(
		double current_vuv, int noise_size,
		const double * const *spectrogram, int fft_size,
		const double * const *aperiodicity, int f0_length, double frame_period,
		double current_time, double fractional_time_shift, int fs,
		const ForwardRealFFT *forward_real_fft, const InverseRealFFT *inverse_real_fft,
		MinimumPhaseAnalysis *minimum_phase, const double *dc_remover,
		double *response, int thread_id = 0
	);

	void getSpectralEnvelope(
		double current_time, double frame_period,
		int f0_length, const double * const *spectrogram, int fft_size,
		double *spectral_envelope
	);

	void getAperiodicRatio(
		double current_time, double frame_period,
		int f0_length, const double * const *aperiodicity, int fft_size,
		double *aperiodic_spectrum
	);

	double getSafeAperiodicity(double x);

	void getPeriodicResponse(
		int fft_size, const double *spectrum,
		const double *aperiodic_ratio, double current_vuv,
		const InverseRealFFT *inverse_real_fft,
		MinimumPhaseAnalysis *minimum_phase, const double *dc_remover,
		double fractional_time_shift, int fs, double *periodic_response
	);

	void getSpectrumWithFractionalTimeShift(
		int fft_size, double coefficient, const InverseRealFFT *inverse_real_fft
	);

	void removeDCComponent(
		const double *periodic_response, int fft_size,
		const double *dc_remover, double *new_periodic_response
	);

	void getAperiodicResponse(
		int noise_size, int fft_size,
		const double *spectrum, const double *aperiodic_ratio, double current_vuv,
		const ForwardRealFFT *forward_real_fft,
		const InverseRealFFT *inverse_real_fft,
		MinimumPhaseAnalysis *minimum_phase, double *aperiodic_response
	);

	void getNoiseSpectrum(
		int noise_size, int fft_size, const ForwardRealFFT *forward_real_fft
	);
};

}

#endif
