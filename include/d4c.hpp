//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise,
// Copyright 2018 Yukara Ikemiya
//
// Band-aperiodicity estimation on the basis of the idea of D4C.
//-----------------------------------------------------------------------------

#ifndef WORLD_CLASS_D4C_HPP
#define WORLD_CLASS_D4C_HPP

#include "world_common.hpp"

namespace world_class
{

typedef struct D4COption{
	double threshold;

	D4COption();
} D4COption;


class D4C
{

public:

	D4C(int fs);
	D4C(int fs, const D4COption &option);
	~D4C();

	void compute(
		const double *x, int x_length,
		const double *temporal_positions, const double *f0, int f0_length,
		int fft_size, double **aperiodicity
	);

private:

	D4COption option_;

	// for D4C
	int fs_;
	int fft_size_d4c_;
	int number_of_aperiodicities_;
	int num_thread_;
	int window_length_;
	double *window_;
	double *coarse_aperiodicity_;
	double *coarse_frequency_axis_;
	ForwardRealFFT* ffts_;

	// for Love Train
	double lowest_f0_;
	int fft_size_lt_;
	ForwardRealFFT* ffts_lt_;

	// input variables
	const double *x_;
	int x_length_;
	const double *time_pos_;
	const double *f0_;
	int f0_length_;
	int fft_size_;

	void prepareForD4c(int fs);
	
	void loveTrain(double *aperiodicity0);
	double loveTrainSub(
		double current_f0, double current_position, int fft_size,
		int boundary0, int boundary1, int boundary2,
		const ForwardRealFFT &forward_real_fft
	);
	void getWindowedWaveform(
		double current_f0, double current_position, int window_type,
		double window_length_ratio, double *waveform
	);

	void generalBody(
		double current_f0, double current_position,
		const ForwardRealFFT &forward_real_fft, double *coarse_aperiodicity
	);

	void getStaticCentroid(
		double current_f0, double current_position,
		const ForwardRealFFT &forward_real_fft, double *static_centroid
	);
	void getCentroid(
		double current_f0, double current_position,
		const ForwardRealFFT &forward_real_fft, double *centroid
	);

	void getSmoothedPowerSpectrum(
		double current_f0, double current_position,
		const ForwardRealFFT &forward_real_fft,
		double *smoothed_power_spectrum
	);

	void getStaticGroupDelay(
		const double *static_centroid, const double *smoothed_power_spectrum,
		double current_f0, double *static_group_delay
	);

	void getCoarseAperiodicity(
		const double *static_group_delay, const ForwardRealFFT &forward_real_fft,
		double *coarse_aperiodicity
	);
    
};

} // end namespace world_class

#endif
