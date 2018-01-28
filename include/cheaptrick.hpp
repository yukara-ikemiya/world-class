//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise,
// Copyright 2018 Yukara Ikemiya
//-----------------------------------------------------------------------------

#ifndef WORLD_CLASS_CHEAPTRICK_HPP
#define WORLD_CLASS_CHEAPTRICK_HPP

#include "world_common.hpp"

namespace world_class
{

typedef struct CheapTrickOption{
	double q1;
	double f0_floor;
	int fft_size;

	CheapTrickOption();
} CheapTrickOption;


class CheapTrick
{

public:

	CheapTrick(int fs);
	CheapTrick(int fs, const CheapTrickOption &option);

	void compute(
		const double *x, int x_length, const double *temporal_positions,
		const double *f0, int f0_length, double **spectrogram
	);

	int getFFTSizeForCheapTrick(int fs, double f0_floor);
  
	double getF0FloorForCheapTrick(int fs, int fft_size);

private:

	CheapTrickOption option_;

	int fs_;
	double f0_floor_;

	const double *x_;
	int x_length_;

	void generalBody(
		double current_f0, double current_position,
		const ForwardRealFFT &forward_real_fft,
		const InverseRealFFT &inverse_real_fft,
		double *spectral_envelope
	);

	void getWindowedWaveform(
		double current_f0, double currnet_position,
		const ForwardRealFFT &forward_real_fft
	);

	void setParametersForGetWindowedWaveform(
		int half_window_length, double currnet_position, double current_f0,
		int *base_index, int *safe_index, double *window
	);

	void getPowerSpectrum(double f0, const ForwardRealFFT &forward_real_fft);

	void addInfinitesimalNoise(const double *input_spectrum, double *output_spectrum);

	void smoothingWithRecovery(
		double f0,
		const ForwardRealFFT &forward_real_fft,
		const InverseRealFFT &inverse_real_fft,
		double *spectral_envelope
	);
};

} // end namespace world_class

#endif
