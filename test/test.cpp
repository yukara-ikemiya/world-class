//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise,
// Copyright 2018 Yukara Ikemiya
//
// test.exe input.wav outout.wav f0 spec
// input.wav  : Input file
//
// output.wav : Output file
// f0         : F0 scaling (a positive number)
// spec       : Formant scaling (a positive number)
//
// Note: This version output three speech synthesized by different algorithms.
//       When the filename is "output.wav", "01output.wav", "02output.wav" and
//       "03output.wav" are generated. They are almost all the same.
//-----------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <chrono>

// For .wav input/output functions.
#include "audioio.hpp"

// WORLD core functions.
#include "world_matlabfunctions.hpp"
#include "harvest.hpp"
#include "cheaptrick.hpp"
#include "d4c.hpp"
#include "synthesis.hpp"

using namespace std;
using namespace world_class;

namespace
{
// time measurement methods
using chrono_tp = chrono::system_clock::time_point;
inline chrono_tp get_time_now() {
	return chrono::system_clock::now();
}

inline double get_elapsed_msec(chrono_tp st, chrono_tp end) {
	return chrono::duration_cast<chrono::nanoseconds>(end-st).count() / 1.0e6;
}
}

//-----------------------------------------------------------------------------
// struct for WORLD
// This struct is an option.
// Users are NOT forced to use this struct.
//-----------------------------------------------------------------------------
typedef struct {
	double frame_period;
	int fs;

	double *f0;
	double *time_axis;
	int f0_length;

	double **spectrogram;
	double **aperiodicity;
	int fft_size;
} WorldParameters;

namespace {

void DisplayInformation(int fs, int nbit, int x_length)
{
	cout << "File information" << endl;
	cout << "Sampling : " << fs << " [Hz] " << nbit << " [Bit]" << endl;
	cout << "Length " << x_length << " [sample]" << endl;
	cout << "Length " << (double)x_length / fs << " [sec]" << endl;
}


void F0EstimationHarvest(
	double *x, int x_length,  WorldParameters *world_parameters)
{
	cout << "\nF0 estimation (Harvest)" << endl;
	
	// You can change the frame period.
	// But the estimation is carried out with 1-ms frame shift.
	HarvestOption option;
	option.frame_period = world_parameters->frame_period;

	// You can set the f0_floor below world::kFloorF0.
	option.f0_floor = 40.0;

	// You can use a Cosine table for fast computation,
	// but the computation is not exactly same as the original algorithm.
	// This option is not included int the original codes.
	
	// option.use_cos_table = true;
	
	auto t1 = get_time_now();
	
	Harvest harvest = Harvest(world_parameters->fs, option);

	auto t2_1 = get_time_now();

	// Parameters setting and memory allocation.
	world_parameters->f0_length = harvest.getSamples(world_parameters->fs, x_length);
	world_parameters->f0 = new double[world_parameters->f0_length];
	world_parameters->time_axis = new double[world_parameters->f0_length];

	auto t2_2 = get_time_now();
	
	harvest.compute(x, x_length, world_parameters->time_axis, world_parameters->f0);

	auto t3 = get_time_now();
	
	cout << "\t initialize:\t" << get_elapsed_msec(t1, t2_1) << " [msec]" << endl;
	cout << "\t compute:\t" << get_elapsed_msec(t2_2, t3) << " [msec]" << endl;
}


void SpectralEnvelopeEstimation(
	double *x, int x_length, WorldParameters *world_parameters)
{
	cout << "\nSpectral envelope estimation (CheapTrick)" << endl;
	
	// Important notice (2017/01/02)
	// You can set the fft_size.
	// Default is GetFFTSizeForCheapTrick(world_parameters->fs, &option);
	// When fft_size changes from default value,
	// a replaced f0_floor will be used in CheapTrick().
	// The lowest F0 that WORLD can work as expected is determined
	// by the following : 3.0 * fs / fft_size.
	CheapTrickOption option;
	option.f0_floor = 71.0;
	// We can directly set fft_size.
	// option.fft_size = 1024;

	// Default value was modified to -0.15.
	// option.q1 = -0.15;
	
	auto t1 = get_time_now();
	
	CheapTrick cheaptrick = CheapTrick(world_parameters->fs, option);

	auto t2_1 = get_time_now();
	
	// Parameters setting and memory allocation.
	int fft_size = cheaptrick.getFFTSizeForCheapTrick(world_parameters->fs, option.f0_floor);
	world_parameters->fft_size = fft_size;
	world_parameters->spectrogram = new double *[world_parameters->f0_length];
	for (int i = 0; i < world_parameters->f0_length; ++i)
		world_parameters->spectrogram[i] =
			new double[world_parameters->fft_size / 2 + 1];
	
	auto t2_2 = get_time_now();
	
	cheaptrick.compute(x, x_length, world_parameters->time_axis,
					   world_parameters->f0, world_parameters->f0_length,
					   world_parameters->spectrogram);

	auto t3 = get_time_now();
	
	cout << "\t initialize:\t" << get_elapsed_msec(t1, t2_1) << " [msec]" << endl;
	cout << "\t compute:\t" << get_elapsed_msec(t2_2, t3) << " [msec]" << endl;
}


void AperiodicityEstimation(
	double *x, int x_length, WorldParameters *world_parameters)
{
	cout << "\nAperiodicity estimation (D4C)" << endl;
	
	// Parameters setting and memory allocation.
	world_parameters->aperiodicity = new double *[world_parameters->f0_length];
	for (int i = 0; i < world_parameters->f0_length; ++i) {
		world_parameters->aperiodicity[i] = new double[world_parameters->fft_size / 2 + 1];
	}
	
	// Comment was modified because it was confusing (2017/12/10).
	// It is used to determine the aperiodicity in whole frequency band.
	// D4C identifies whether the frame is voiced segment even if it had an F0.
	// If the estimated value falls below the threshold,
	// the aperiodicity in whole frequency band will set to 1.0.
	// If you want to use the conventional D4C, please set the threshold to 0.0.
	D4COption option;
	option.threshold = 0.85;
	
	auto t1 = get_time_now();
	
	D4C d4c = D4C(world_parameters->fs, option);

	auto t2 = get_time_now();
	
	d4c.compute(x, x_length, world_parameters->time_axis,
				world_parameters->f0, world_parameters->f0_length,
				world_parameters->fft_size, world_parameters->aperiodicity);
		
	auto t3 = get_time_now();
	
	cout << "\t initialize:\t" << get_elapsed_msec(t1, t2) << " [msec]" << endl;
	cout << "\t compute:\t" << get_elapsed_msec(t2, t3) << " [msec]" << endl;
}


void ParameterModification(
	int argc, char *argv[], int fs, int f0_length,
	int fft_size, double *f0, double **spectrogram)
{
	// F0 scaling
	if (argc >= 4) {
		double shift = atof(argv[3]);
		for (int i = 0; i < f0_length; ++i) f0[i] *= shift;
	}
	
	if (argc < 5) return;

	// Spectral stretching
	double ratio = atof(argv[4]);
	double *freq_axis1 = new double[fft_size];
	double *freq_axis2 = new double[fft_size];
	double *spectrum1 = new double[fft_size];
	double *spectrum2 = new double[fft_size];

	for (int i = 0; i <= fft_size / 2; ++i) {
		freq_axis1[i] = ratio * i / fft_size * fs;
		freq_axis2[i] = static_cast<double>(i) / fft_size * fs;
	}
	
	for (int i = 0; i < f0_length; ++i) {
		for (int j = 0; j <= fft_size / 2; ++j)
			spectrum1[j] = log(spectrogram[i][j]);
		interp1(freq_axis1, spectrum1, fft_size / 2 + 1, freq_axis2,
				fft_size / 2 + 1, spectrum2);
		for (int j = 0; j <= fft_size / 2; ++j)
			spectrogram[i][j] = exp(spectrum2[j]);
		if (ratio >= 1.0) continue;
		for (int j = static_cast<int>(fft_size / 2.0 * ratio);
			 j <= fft_size / 2; ++j)
			spectrogram[i][j] =
				spectrogram[i][static_cast<int>(fft_size / 2.0 * ratio) - 1];
	}
	
	delete[] spectrum1;
	delete[] spectrum2;
	delete[] freq_axis1;
	delete[] freq_axis2;
}

void WaveformSynthesis1(
	WorldParameters *world_parameters, int fs, int y_length, double *y)
{
	cout << "\nSynthesis 1 (conventional algorithm)" << endl;
	
	auto t1 = get_time_now();
	
	Synthesis synthesis = Synthesis(fs, world_parameters->fft_size, world_parameters->frame_period);
	
	auto t2 = get_time_now();
	
	synthesis.compute(world_parameters->f0, world_parameters->f0_length,
					  world_parameters->spectrogram, world_parameters->aperiodicity,
					  y_length, y);
	
	auto t3 = get_time_now();
	
	cout << "\t initialize:\t" << get_elapsed_msec(t1, t2) << " [msec]" << endl;
	cout << "\t compute:\t" << get_elapsed_msec(t2, t3) << " [msec]" << endl;
}

void DestroyMemory(WorldParameters *world_parameters) {
	delete[] world_parameters->time_axis;
	delete[] world_parameters->f0;
	for (int i = 0; i < world_parameters->f0_length; ++i) {
		delete[] world_parameters->spectrogram[i];
		delete[] world_parameters->aperiodicity[i];
	}
	delete[] world_parameters->spectrogram;
	delete[] world_parameters->aperiodicity;
}

}  // end unnamed namespace


//-----------------------------------------------------------------------------
// Test program.
// test.exe input.wav outout f0 spec flag
// input.wav  : argv[1] Input file path
// output	  : argv[2] Output file name
// f0         : argv[3] F0 scaling (a positive number)
// spec       : argv[4] Formant shift (a positive number)
//-----------------------------------------------------------------------------
int main(int argc, char *argv[]) {
	if (argc != 2 && argc != 3 && argc != 4 && argc != 5) {
		cout << "Usage : ./test [input wav path] [output wav name]" << endl;
		cout << "Example: ./test ./aiueo.wav output" << endl;
		return -2;
	}

	// Memory allocation is carried out in advanse.
	// This is for compatibility with C language.
	int x_length = GetAudioLength(argv[1]);
	if (x_length <= 0) {
		if (x_length == 0) { cerr << "error: File not found.\n" << endl; }
		else { cerr << "error: The file is not .wav format.\n" << endl; }
		return -1;
	}
	
	
	// wavread() must be called after GetAudioLength().
	int fs, nbit;
	double *x = new double[x_length];
	wavread(argv[1], &fs, &nbit, x);
	
	DisplayInformation(fs, nbit, x_length);

	//---------------------------------------------------------------------------
	// Analysis part
	//---------------------------------------------------------------------------
	WorldParameters world_parameters;
	// You must set fs and frame_period before analysis/synthesis.
	world_parameters.fs = fs;
	// 5.0 ms is the default value.
	world_parameters.frame_period = 5.0;

	// F0 estimation

	// Harvest
	F0EstimationHarvest(x, x_length, &world_parameters);

	double f0_sum = 0;
	for (int i = 0; i < world_parameters.f0_length; i++) {
		f0_sum += world_parameters.f0[i];
	}
    
	// Spectral envelope estimation
	SpectralEnvelopeEstimation(x, x_length, &world_parameters);

	double spec_sum = 0;
	for (int i = 0; i < world_parameters.f0_length; ++i)
		for (int j = 0; j < world_parameters.fft_size / 2 + 1; ++j)
		{ spec_sum +=  world_parameters.spectrogram[i][j]; }
  
	// Aperiodicity estimation by D4C
	AperiodicityEstimation(x, x_length, &world_parameters);

	double ap_sum = 0;
	for (int i = 0; i < world_parameters.f0_length; ++i)
		for (int j = 0; j < world_parameters.fft_size / 2 + 1; ++j)
		{ ap_sum +=  world_parameters.aperiodicity[i][j]; }
	
	// Note that F0 must not be changed until all parameters are estimated.
	ParameterModification(argc, argv, fs, world_parameters.f0_length,
						  world_parameters.fft_size, world_parameters.f0,
						  world_parameters.spectrogram);

	
	//---------------------------------------------------------------------------
	// Synthesis part
	// There are three samples in speech synthesis
	// 1: Conventional synthesis
	// 2: Example of real-time synthesis
	// 3: Example of real-time synthesis (Ring buffer is efficiently used)
	//---------------------------------------------------------------------------
	char filename[100];
	// The length of the output waveform
	int y_length = static_cast<int>((world_parameters.f0_length - 1) *
									world_parameters.frame_period / 1000.0 * fs) + 1;
	double *y = new double[y_length]();

	
	// Synthesis 1 (conventional synthesis)
	WaveformSynthesis1(&world_parameters, fs, y_length, y);
	if (argc == 2) { sprintf(filename, "output_1.wav", argv[2]); }
	else { sprintf(filename, "%s_1.wav", argv[2]); }
	wavwrite(y, y_length, fs, 16, filename);

	double synt_sum_1 = 0;
	for (int ii = 0; ii < y_length; ii++) {
		synt_sum_1 += y[ii];
	}

	delete[] y;
	delete[] x;
	DestroyMemory(&world_parameters);

	cout << "complete." << endl;
	return 0;
}
