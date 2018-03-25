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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <iostream>

#if (defined (__WIN32__) || defined (_WIN32)) && !defined (__MINGW32__)
#include <conio.h>
#include <windows.h>
#pragma comment(lib, "winmm.lib")
#pragma warning(disable : 4996)
#endif
#if (defined (__linux__) || defined(__CYGWIN__) || defined(__APPLE__))
#include <stdint.h>
#include <sys/time.h>
#endif

// For .wav input/output functions.
#include "audioio.hpp"

// WORLD core functions.
#include "world_matlabfunctions.hpp"
#include "harvest.hpp"
#include "cheaptrick.hpp"
#include "d4c.hpp"
// #include "world/stonemask.h"
#include "synthesis.hpp"
// #include "world/synthesisrealtime.h"

#if (defined (__linux__) || defined(__CYGWIN__) || defined(__APPLE__))
// Linux porting section: implement timeGetTime() by gettimeofday(),
#ifndef DWORD
#define DWORD uint32_t
#endif
DWORD timeGetTime() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	DWORD ret = static_cast<DWORD>(tv.tv_usec / 1000 + tv.tv_sec * 1000);
	return ret;
}
#endif

using namespace std;
using namespace world_class;

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
	printf("File information\n");
	printf("Sampling : %d Hz %d Bit\n", fs, nbit);
	printf("Length %d [sample]\n", x_length);
	printf("Length %f [sec]\n", static_cast<double>(x_length) / fs);
}


void F0EstimationHarvest(double *x, int x_length,
						 WorldParameters *world_parameters)
{
	HarvestOption option;

	// You can change the frame period.
	// But the estimation is carried out with 1-ms frame shift.
	option.frame_period = world_parameters->frame_period;

	// You can set the f0_floor below world::kFloorF0.
	option.f0_floor = 40.0;

	// Harvest class
	Harvest harvest = Harvest(option);

	// Parameters setting and memory allocation.
	world_parameters->f0_length = harvest.getSamples(world_parameters->fs, x_length);
	world_parameters->f0 = new double[world_parameters->f0_length];
	world_parameters->time_axis = new double[world_parameters->f0_length];

	printf("\nAnalysis\n");
	DWORD elapsed_time = timeGetTime();
	harvest.compute(x, x_length, world_parameters->fs,
					world_parameters->time_axis, world_parameters->f0);
	printf("Harvest: %d [msec]\n", timeGetTime() - elapsed_time);
}


void SpectralEnvelopeEstimation(double *x, int x_length,
								WorldParameters *world_parameters)
{
	CheapTrickOption option;
  
	// Default value was modified to -0.15.
	// option.q1 = -0.15;

	// Important notice (2017/01/02)
	// You can set the fft_size.
	// Default is GetFFTSizeForCheapTrick(world_parameters->fs, &option);
	// When fft_size changes from default value,
	// a replaced f0_floor will be used in CheapTrick().
	// The lowest F0 that WORLD can work as expected is determined
	// by the following : 3.0 * fs / fft_size.
	option.f0_floor = 71.0;
	// We can directly set fft_size.
	// option.fft_size = 1024;

	// CheapTrick class
	CheapTrick cheaptrick = CheapTrick(world_parameters->fs, option);
  
	int fft_size = cheaptrick.getFFTSizeForCheapTrick(world_parameters->fs, option.f0_floor);
	// Parameters setting and memory allocation.
	world_parameters->fft_size = fft_size;
	world_parameters->spectrogram = new double *[world_parameters->f0_length];
	for (int i = 0; i < world_parameters->f0_length; ++i)
		world_parameters->spectrogram[i] =
			new double[world_parameters->fft_size / 2 + 1];

	DWORD elapsed_time = timeGetTime();
	cheaptrick.compute(x, x_length, world_parameters->time_axis,
					   world_parameters->f0, world_parameters->f0_length,
					   world_parameters->spectrogram);
	printf("CheapTrick: %d [msec]\n", timeGetTime() - elapsed_time);
}


void AperiodicityEstimation(double *x, int x_length,
							WorldParameters *world_parameters)
{
	D4COption option;

	// Comment was modified because it was confusing (2017/12/10).
	// It is used to determine the aperiodicity in whole frequency band.
	// D4C identifies whether the frame is voiced segment even if it had an F0.
	// If the estimated value falls below the threshold,
	// the aperiodicity in whole frequency band will set to 1.0.
	// If you want to use the conventional D4C, please set the threshold to 0.0.
	option.threshold = 0.85;

	// D4C class
	D4C d4c = D4C(world_parameters->fs, option);

	// Parameters setting and memory allocation.
	world_parameters->aperiodicity = new double *[world_parameters->f0_length];
	for (int i = 0; i < world_parameters->f0_length; ++i) {
		world_parameters->aperiodicity[i] = new double[world_parameters->fft_size / 2 + 1];
	}

	DWORD elapsed_time = timeGetTime();
	d4c.compute(x, x_length, world_parameters->time_axis,
				world_parameters->f0, world_parameters->f0_length,
				world_parameters->fft_size, world_parameters->aperiodicity);
	printf("D4C: %d [msec]\n", timeGetTime() - elapsed_time);
}


void ParameterModification(int argc, char *argv[], int fs, int f0_length,
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

void WaveformSynthesis1(WorldParameters *world_parameters, int fs,
						int y_length, double *y)
{
	DWORD elapsed_time;
	// Synthesis by the aperiodicity
	printf("\nSynthesis 1 (conventional algorithm)\n");
	elapsed_time = timeGetTime();
	Synthesis synthesis = Synthesis(fs, world_parameters->fft_size, world_parameters->frame_period);
	synthesis.compute(world_parameters->f0, world_parameters->f0_length,
					  world_parameters->spectrogram, world_parameters->aperiodicity,
					  y_length, y);
	printf("WORLD: %d [msec]\n", timeGetTime() - elapsed_time);
}

/*
void WaveformSynthesis2(WorldParameters *world_parameters, int fs,
						int y_length, double *y) {
	DWORD elapsed_time;
	printf("\nSynthesis 2 (All frames are added at the same time)\n");
	elapsed_time = timeGetTime();

	WorldSynthesizer synthesizer = { 0 };
	int buffer_size = 64;
	InitializeSynthesizer(world_parameters->fs, world_parameters->frame_period,
						  world_parameters->fft_size, buffer_size, 1, &synthesizer);

	// All parameters are added at the same time.
	AddParameters(world_parameters->f0, world_parameters->f0_length,
				  world_parameters->spectrogram, world_parameters->aperiodicity,
				  &synthesizer);

	int index;
	for (int i = 0; Synthesis2(&synthesizer) != 0; ++i) {
		index = i * buffer_size;
		for (int j = 0; j < buffer_size; ++j)
			y[j + index] = synthesizer.buffer[j];
	}

	printf("WORLD: %d [msec]\n", timeGetTime() - elapsed_time);
	DestroySynthesizer(&synthesizer);
}
*/

  /*
  void WaveformSynthesis3(WorldParameters *world_parameters, int fs,
  int y_length, double *y) {
  DWORD elapsed_time;
  // Synthesis by the aperiodicity
  printf("\nSynthesis 3 (Ring buffer is efficiently used.)\n");
  elapsed_time = timeGetTime();

  WorldSynthesizer synthesizer = { 0 };
  int buffer_size = 64;
  InitializeSynthesizer(world_parameters->fs, world_parameters->frame_period,
  world_parameters->fft_size, buffer_size, 100, &synthesizer);

  int offset = 0;
  int index = 0;
  for (int i = 0; i < world_parameters->f0_length;) {
  // Add one frame (i shows the frame index that should be added)
  if (AddParameters(&world_parameters->f0[i], 1,
  &world_parameters->spectrogram[i], &world_parameters->aperiodicity[i],
  &synthesizer) == 1) ++i;

  // Synthesize speech with length of buffer_size sample.
  // It is repeated until the function returns 0
  // (it suggests that the synthesizer cannot generate speech).
  while (Synthesis2(&synthesizer) != 0) {
  index = offset * buffer_size;
  for (int j = 0; j < buffer_size; ++j)
  y[j + index] = synthesizer.buffer[j];
  offset++;
  }

  // Check the "Lock" (Please see synthesisrealtime.h)
  if (IsLocked(&synthesizer) == 1) {
  printf("Locked!\n");
  break;
  }
  }

  printf("WORLD: %d [msec]\n", timeGetTime() - elapsed_time);
  DestroySynthesizer(&synthesizer);
  }
*/

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

}  // namespace

//-----------------------------------------------------------------------------
// Test program.
// test.exe input.wav outout.wav f0 spec flag
// input.wav  : argv[1] Input file
// output.wav : argv[2] Output file
// f0         : argv[3] F0 scaling (a positive number)
// spec       : argv[4] Formant shift (a positive number)
//-----------------------------------------------------------------------------
int main(int argc, char *argv[]) {
	if (argc != 2 && argc != 3 && argc != 4 && argc != 5) {
		printf("error\n");
		return -2;
	}

	// Memory allocation is carried out in advanse.
	// This is for compatibility with C language.
	int x_length = GetAudioLength(argv[1]);
	if (x_length <= 0) {
		if (x_length == 0) printf("error: File not found.\n");
		else printf("error: The file is not .wav format.\n");
		return -1;
	}
	double *x = new double[x_length];
	// wavread() must be called after GetAudioLength().
	int fs, nbit;
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
	// DIO
	// F0EstimationDio(x, x_length, &world_parameters);

	// Harvest
	F0EstimationHarvest(x, x_length, &world_parameters);

	double f0_sum = 0;
	for (int i = 0; i < world_parameters.f0_length; i++) {
		f0_sum += world_parameters.f0[i];
	}

	cout << "f0 sum (Harvest): " << f0_sum << endl;
    
	// Spectral envelope estimation
	SpectralEnvelopeEstimation(x, x_length, &world_parameters);

	double spec_sum = 0;
	for (int i = 0; i < world_parameters.f0_length; ++i)
		for (int j = 0; j < world_parameters.fft_size / 2 + 1; ++j)
		{ spec_sum +=  world_parameters.spectrogram[i][j]; }

	cout << "spec sum (CheapTrick): " << spec_sum << endl;
  
	// Aperiodicity estimation by D4C
	AperiodicityEstimation(x, x_length, &world_parameters);

	double ap_sum = 0;
	for (int i = 0; i < world_parameters.f0_length; ++i)
		for (int j = 0; j < world_parameters.fft_size / 2 + 1; ++j)
		{ ap_sum +=  world_parameters.aperiodicity[i][j]; }

	cout << "ap_sum (D4C): " << ap_sum << endl;
	
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
	sprintf(filename, "01%s", argv[2]);
	wavwrite(y, y_length, fs, 16, filename);

	double synt_sum_1 = 0;
	for (int ii = 0; ii < y_length; ii++) {
		synt_sum_1 += y[ii];
	}

	cout << "synthesis1_sum: " << synt_sum_1 << endl;

	/*
	// Synthesis 2 (All frames are added at the same time)
	for (int i = 0; i < y_length; ++i) y[i] = 0.0;
	WaveformSynthesis2(&world_parameters, fs, y_length, y);
	sprintf(filename, "02%s", argv[2]);
	wavwrite(y, y_length, fs, 16, filename);

	// Synthesis 3 (Ring buffer is efficiently used.)
	for (int i = 0; i < y_length; ++i) y[i] = 0.0;
	WaveformSynthesis3(&world_parameters, fs, y_length, y);
	sprintf(filename, "03%s", argv[2]);
	wavwrite(y, y_length, fs, 16, filename);
	*/

	// delete[] y;
	delete[] x;
	DestroyMemory(&world_parameters);

	printf("complete.\n");
	return 0;
}
