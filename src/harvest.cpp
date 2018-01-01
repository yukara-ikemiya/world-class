#include "harvest.hpp"

#include <cmath>
#include <algorithm>
#include <numeric>
#include <cstring>
#include <iostream>

#include "world_common.h"
#include "world_constantnumbers.h"
#include "world_fft.h"
#include "world_matlabfunctions.h"


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

#if (defined (__linux__) || defined(__CYGWIN__) || defined(__APPLE__))
// Linux porting section: implement timeGetTime() by gettimeofday(),
#ifndef DWORD
#define DWORD uint32_t
#endif
static DWORD timeGetTime() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  DWORD ret = static_cast<DWORD>(tv.tv_usec / 1000 + tv.tv_sec * 1000);
  return ret;
}
#endif

using namespace std;


namespace world_class
{

HarvestOption::HarvestOption()
  : f0_floor(world::kFloorF0), f0_ceil(world::kCeilF0), frame_period(5)
{}


Harvest::Harvest()
{}

Harvest::Harvest(const HarvestOption& option)
{
  option_.f0_floor = option.f0_floor;
  option_.f0_ceil = option.f0_ceil;
  option_.frame_period = option.frame_period;
}

int Harvest::getSamples(const int fs, const int x_length, const double frame_period)
{
  return static_cast<int>(1000.0 * x_length / fs / frame_period) + 1;
}

int Harvest::getSamples(const int fs, const int x_length)
{
  return static_cast<int>(1000.0 * x_length / fs / option_.frame_period) + 1;
}

void Harvest::compute(const double* x, const int x_length, const int fs,
		      double *temporal_positions, double *f0)
{
  // Several parameters will be controllable for debug.
  double target_fs = 8000.0;
  int dimension_ratio = matlab_round(fs / target_fs);
  double channels_in_octave = 40;

  if (option_.frame_period == 1.0) {
    generalBody(x, x_length, fs, 1,
		channels_in_octave, dimension_ratio,
		temporal_positions, f0);
    return;
  }

  int basic_frame_period = 1;
  int basic_f0_length = getSamples(fs, x_length, basic_frame_period);
  double *basic_f0 = new double[basic_f0_length];
  double *basic_temporal_positions = new double[basic_f0_length];
  
  generalBody(x, x_length, fs, basic_frame_period,
		     channels_in_octave, dimension_ratio,
		     basic_temporal_positions, basic_f0);

  int f0_length = getSamples(fs, x_length, option_.frame_period);
  for (int i = 0; i < f0_length; ++i) {
    temporal_positions[i] = i * option_.frame_period / 1000.0;
    f0[i] = basic_f0[MyMinInt(basic_f0_length - 1,
			      matlab_round(temporal_positions[i] * 1000.0))];
  }

  delete[] basic_f0;
  delete[] basic_temporal_positions;
}


//-----------------------------------------------------------------------------
// Since the waveform of beginning and ending after decimate include noise,
// the input waveform is extended. This is the processing for the
// compatibility with MATLAB version.
//-----------------------------------------------------------------------------
void Harvest::getWaveformAndSpectrumSub(const double *x, const int x_length,
					const int y_length, const double actual_fs,
					const int decimation_ratio,
					double *y)
{
  if (decimation_ratio == 1) {
    copy(x, x + x_length, y);
    return;
  }

  int lag = static_cast<int>(ceil(140.0 / decimation_ratio) * decimation_ratio);
  int new_x_length = x_length + lag * 2;
  double *new_y = new double[new_x_length]();
  double *new_x = new double[new_x_length];
  
  for (int i = 0; i < lag; ++i) { new_x[i] = x[0]; }
  copy(x, x + x_length, new_x + lag);
  for (int i = lag + x_length; i < new_x_length; ++i) { new_x[i] = x[x_length - 1]; }

  decimate(new_x, new_x_length, decimation_ratio, new_y);
  for (int i = 0; i < y_length; ++i) { y[i] = new_y[lag / decimation_ratio + i]; }

  delete[] new_x;
  delete[] new_y;
}

//-----------------------------------------------------------------------------
// GetWaveformAndSpectrum() calculates the downsampled signal and its spectrum
//-----------------------------------------------------------------------------
void Harvest::getWaveformAndSpectrum(const double *x, const int x_length,
				     const int y_length, const double actual_fs,
				     const int fft_size, const int decimation_ratio,
				     double *y, fft_complex *y_spectrum)
{
  
  // Initialization
  // memset(y, 0, sizeof(double) * fft_size); // initialized before this function

  // Processing for the compatibility with MATLAB version
  getWaveformAndSpectrumSub(x, x_length, y_length, actual_fs,
			    decimation_ratio, y);

  // Removal of the DC component (y = y - mean value of y)
  // double mean_y = 0.0;
  // for (int i = 0; i < y_length; ++i) mean_y += y[i];
  double mean_y = accumulate(y, y + y_length, 0);
  mean_y /= y_length;
  for (int i = 0; i < y_length; ++i) y[i] -= mean_y;
  for (int i = y_length; i < fft_size; ++i) y[i] = 0.0;

  fft_plan forwardFFT = fft_plan_dft_r2c_1d(fft_size, y, y_spectrum, FFT_ESTIMATE);
  fft_execute(forwardFFT);

  fft_destroy_plan(forwardFFT);
}

//-----------------------------------------------------------------------------
// DestroyZeroCrossings() frees the memory of array in the struct
//-----------------------------------------------------------------------------
// static void DestroyZeroCrossings(ZeroCrossings *zero_crossings) {
//   delete[] zero_crossings.negative_interval_locations;
//   delete[] zero_crossings.positive_interval_locations;
//   delete[] zero_crossings.peak_interval_locations;
//   delete[] zero_crossings.dip_interval_locations;
//   delete[] zero_crossings.negative_intervals;
//   delete[] zero_crossings.positive_intervals;
//   delete[] zero_crossings.peak_intervals;
//   delete[] zero_crossings.dip_intervals;
// }


//-----------------------------------------------------------------------------
// SearchF0Base() gets the F0 with the highest score.
//-----------------------------------------------------------------------------
void Harvest::searchF0Base(const double * const *f0_candidates, const double * const *f0_scores,
			   const int f0_length, const int number_of_candidates,
			   double *base_f0_contour)
{
  double tmp_best_score;
  
  for (int i = 0; i < f0_length; ++i) {
    const double* ptr_f0_scores = f0_scores[i];
    const double* ptr_f0_candidates = f0_candidates[i];
    
    base_f0_contour[i] = tmp_best_score = 0.0;
    for (int j = 0; j < number_of_candidates; ++j) {
      if (ptr_f0_scores[j] > tmp_best_score) {
  	base_f0_contour[i] = ptr_f0_candidates[j];
  	tmp_best_score = ptr_f0_scores[j];
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Step 1: Rapid change of F0 contour is replaced by 0.
//-----------------------------------------------------------------------------
void Harvest::fixStep1(const double *f0_base, const int f0_length,
		       const double allowed_range, double *f0_step1)
{
  f0_step1[0] = f0_step1[1] = 0.0;
  double reference_f0;
  
  for (int i = 2; i < f0_length; ++i) {
    if (f0_base[i] == 0.0) { continue; }
    
    reference_f0 = f0_base[i - 1] * 2 - f0_base[i - 2];
    f0_step1[i] =
      (fabs((f0_base[i] - reference_f0) / reference_f0) > allowed_range &&
       fabs((f0_base[i] - f0_base[i - 1])) / f0_base[i - 1] > allowed_range)
      ? 0.0 : f0_base[i];
  }
}

//-----------------------------------------------------------------------------
// GetBoundaryList() detects boundaries between voiced and unvoiced sections.
//-----------------------------------------------------------------------------
int Harvest::getBoundaryList(const double *f0, const int f0_length, int *boundary_list)
{
  int number_of_boundaries = 0;
  int *vuv = new int[f0_length];

  vuv[0] = vuv[f0_length - 1] = 0;
  for (int i = 1; i < f0_length - 1; ++i) { vuv[i] = f0[i] > 0 ? 1 : 0; }
  
  for (int i = 1; i < f0_length; ++i) {
    if (vuv[i] - vuv[i - 1] != 0) {
      boundary_list[number_of_boundaries] = i - number_of_boundaries % 2;
      number_of_boundaries++;
    }
  }

  delete[] vuv;
  
  return number_of_boundaries;
}

//-----------------------------------------------------------------------------
// Step 2: Voiced sections with a short period are removed.
//-----------------------------------------------------------------------------
void Harvest::fixStep2(const double *f0_step1, const int f0_length,
		       const int voice_range_minimum, double *f0_step2)
{
  copy(f0_step1, f0_step1 + f0_length, f0_step2);

  int *boundary_list = new int[f0_length];
  int number_of_boundaries = getBoundaryList(f0_step1, f0_length, boundary_list);

  for (int i = 0; i < number_of_boundaries / 2; ++i) {
    if (boundary_list[i * 2 + 1] - boundary_list[i * 2] >= voice_range_minimum)
      { continue; }
    for (int j = boundary_list[i * 2]; j <= boundary_list[(i * 2) + 1]; ++j)
      { f0_step2[j] = 0.0; }
  }
  
  delete[] boundary_list;
}

//-----------------------------------------------------------------------------
// abs() often causes bugs, an original function is used.
//-----------------------------------------------------------------------------
inline int Harvest::MyAbsInt(int x)
{
  return x > 0 ? x : -x;
}

//-----------------------------------------------------------------------------
// SelectBestF0() obtains the nearlest F0 in reference_f0.
//-----------------------------------------------------------------------------
double Harvest::selectBestF0(const double reference_f0, const double *f0_candidates,
			     const int number_of_candidates, const double allowed_range,
			     double &best_error)
{
  double best_f0 = 0.0;
  best_error = allowed_range;

  double tmp;
  for (int i = 0; i < number_of_candidates; ++i) {
    tmp = fabs(reference_f0 - f0_candidates[i]) / reference_f0;
    
    if (tmp > best_error) continue;
    
    best_f0 = f0_candidates[i];
    best_error = tmp;
  }

  return best_f0;
}

//-----------------------------------------------------------------------------
// ExtendF0() : The Hand erasing the Space.
// The subfunction of Extend().
//-----------------------------------------------------------------------------
int Harvest::extendF0(const double *f0, const int f0_length, const int origin,
		      const int last_point, const int shift, const double * const *f0_candidates,
		      const int number_of_candidates, const double allowed_range, double *extended_f0)
{
  int threshold = 4;
  double tmp_f0 = extended_f0[origin];
  int shifted_origin = origin;

  int distance = MyAbsInt(last_point - origin);
  int *index_list = new int[distance + 1];
  for (int i = 0; i <= distance; ++i) { index_list[i] = origin + shift * i; }

  int count = 0;
  double dammy;
  for (int i = 0; i <= distance; ++i) {
    extended_f0[index_list[i] + shift] = selectBestF0(tmp_f0, f0_candidates[index_list[i] + shift],
						      number_of_candidates, allowed_range, dammy);
    
    if (extended_f0[index_list[i] + shift] == 0.0) {
      count++;
    } else {
      tmp_f0 = extended_f0[index_list[i] + shift];
      count = 0;
      shifted_origin = index_list[i] + shift;
    }
    
    if (count == threshold) break;
  }

  delete[] index_list;
  
  return shifted_origin;
}

//-----------------------------------------------------------------------------
// Swap the f0 contour and boundary.
// It is used in ExtendSub() and MergeF0();
//-----------------------------------------------------------------------------
inline void Harvest::swapArray(const int index1, const int index2, double **f0, int *boundary)
{
  double *tmp_pointer;
  int tmp_index;
  tmp_pointer = f0[index1];
  f0[index1] = f0[index2];
  f0[index2] = tmp_pointer;
  tmp_index = boundary[index1 * 2];
  boundary[index1 * 2] = boundary[index2 * 2];
  boundary[index2 * 2] = tmp_index;
  tmp_index = boundary[index1 * 2 + 1];
  boundary[index1 * 2 + 1] = boundary[index2 * 2 + 1];
  boundary[index2 * 2 + 1] = tmp_index;
}

// int Harvest::extendSub(const double * const *extended_f0,
// 		     const int *boundary_list, int number_of_sections,
// 		     double **selected_extended_f0, int *selected_boundary_list) {
//   double threshold = 2200.0;
//   int count = 0;
//   double mean_f0 = 0.0;
//   int st, ed;
//   for (int i = 0; i < number_of_sections; ++i) {
//     st = boundary_list[i * 2];
//     ed = boundary_list[i * 2 + 1];
//     for (int j = st; j < ed; ++j) mean_f0 += extended_f0[i][j];
//     mean_f0 /= ed - st;
//     if (threshold / mean_f0 < ed - st)
//       swapArray(count++, i, selected_extended_f0, selected_boundary_list);
//   }
//   return count;
// }

//-----------------------------------------------------------------------------
// Extend() : The Hand erasing the Space.
//-----------------------------------------------------------------------------
int Harvest::extend(const double * const *multi_channel_f0,
		    const int number_of_sections, const int f0_length, const int *boundary_list,
		    const double * const *f0_candidates, const int number_of_candidates,
		    const double allowed_range, double **extended_f0, int *shifted_boundary_list)
{
  int threshold = 100;
  for (int i = 0; i < number_of_sections; ++i) {
    shifted_boundary_list[i * 2 + 1] =
      extendF0(multi_channel_f0[i],
	       f0_length, boundary_list[i * 2 + 1],
	       MyMinInt(f0_length - 2, boundary_list[i * 2 + 1] + threshold), 1,
	       f0_candidates, number_of_candidates, allowed_range, extended_f0[i]);
    shifted_boundary_list[i * 2] =
      extendF0(multi_channel_f0[i], f0_length,
	       boundary_list[i * 2], MyMaxInt(1, boundary_list[i * 2] - threshold), -1,
	       f0_candidates, number_of_candidates, allowed_range, extended_f0[i]);
  }

  // extendSub
  double threshold2 = 2200.0;
  int count = 0;
  double mean_f0 = 0.0;
  int st, ed;
  for (int i = 0; i < number_of_sections; ++i) {
    st = shifted_boundary_list[i * 2];
    ed = shifted_boundary_list[i * 2 + 1];
    for (int j = st; j < ed; ++j) mean_f0 += multi_channel_f0[i][j];
    mean_f0 /= ed - st;
    if (threshold2 / mean_f0 < ed - st)
      swapArray(count++, i, extended_f0, shifted_boundary_list);
  }

  // return extendSub(multi_channel_f0, shifted_boundary_list,
  // 		   number_of_sections, extended_f0, shifted_boundary_list);
  return count;
}

//-----------------------------------------------------------------------------
// Indices are sorted.
//-----------------------------------------------------------------------------
// static void MakeSortedOrder(const int *boundary_list, int number_of_sections,
// 			    int *order) {
//   for (int i = 0; i < number_of_sections; ++i) order[i] = i;
//   int tmp;
//   for (int i = 1; i < number_of_sections; ++i)
//     for (int j = i - 1; j >= 0; --j)
//       if (boundary_list[order[j] * 2] > boundary_list[order[i] * 2]) {
// 	tmp = order[i];
// 	order[i] = order[j];
// 	order[j] = tmp;
//       } else {
// 	break;
//       }
// }

//-----------------------------------------------------------------------------
// Serach the highest score with the candidate F0.
//-----------------------------------------------------------------------------
inline double Harvest::searchScore(const double f0, const double *f0_candidates,
				   const double  *f0_scores, const int number_of_candidates)
{
  double score = 0.0;
  for (int i = 0; i < number_of_candidates; ++i)
    if (f0 == f0_candidates[i] && score < f0_scores[i]) { score = f0_scores[i]; }
  return score;
}

//-----------------------------------------------------------------------------
// Subfunction of MergeF0()
//-----------------------------------------------------------------------------
int Harvest::mergeF0Sub(const double *f0_1, const int f0_length, const int st1, const int ed1,
			const double *f0_2, const int st2, const int ed2,
			const double * const *f0_candidates,
			const double * const *f0_scores, const int number_of_candidates,
			double *merged_f0)
{
  if (st1 <= st2 && ed1 >= ed2) return ed1;

  double score1 = 0.0;
  double score2 = 0.0;
  for (int i = st2; i <= ed1; ++i) {
    score1 += searchScore(f0_1[i], f0_candidates[i], f0_scores[i],
			  number_of_candidates);
    score2 += searchScore(f0_2[i], f0_candidates[i], f0_scores[i],
			  number_of_candidates);
  }
  
  if (score1 > score2)
    { copy(f0_2 + ed1, f0_2 + ed2 + 1, merged_f0 + ed1); }
  else
    { copy(f0_2 + st2, f0_2 + ed2 + 1, merged_f0 + st2); }

  return ed2;
}

//-----------------------------------------------------------------------------
// Overlapped F0 contours are merged by the likability score.
//-----------------------------------------------------------------------------
void Harvest::mergeF0(const double * const *multi_channel_f0, int *boundary_list,
		      const int number_of_channels, const int f0_length,
		      const double * const *f0_candidates,
		      const double * const *f0_scores, const int number_of_candidates,
		      double *merged_f0)
{
  // argsort
  int *order = new int[number_of_channels];
  int n = 0;
  generate(order, order + number_of_channels, [&]{ return n++; });
  sort(order, order + number_of_channels,
       [&](int i1, int i2) { return boundary_list[i1 * 2] <= boundary_list[i2 * 2]; } );
  
  // MakeSortedOrder(boundary_list, number_of_channels, order);

  copy(multi_channel_f0[0], multi_channel_f0[0] + f0_length, merged_f0);

  for (int i = 1; i < number_of_channels; ++i) {
    if (boundary_list[order[i] * 2] - boundary_list[1] > 0) {
      // for (int j = boundary_list[order[i] * 2];
      // 	   j <= boundary_list[order[i] * 2 + 1]; ++j)
      // 	merged_f0[j] = multi_channel_f0[order[i]][j];
      int index1 = boundary_list[order[i] * 2];
      int index2 = boundary_list[order[i] * 2 + 1];
      copy(multi_channel_f0[order[i]] + index1,
	   multi_channel_f0[order[i]] + index2 + 1,
	   merged_f0);
      boundary_list[0] = index1;
      boundary_list[1] = index2;
    } else {
      boundary_list[1] =
	mergeF0Sub(merged_f0, f0_length, boundary_list[0], boundary_list[1],
		   multi_channel_f0[order[i]], boundary_list[order[i] * 2],
		   boundary_list[order[i] * 2 + 1], f0_candidates, f0_scores,
		   number_of_candidates, merged_f0);
    }
  }

  delete[] order;
}


//-----------------------------------------------------------------------------
// GetMultiChannelF0() separates each voiced section into independent channel.
//-----------------------------------------------------------------------------
void Harvest::getMultiChannelF0(const double *f0, const int f0_length,
				const int *boundary_list, const int number_of_boundaries,
				double **multi_channel_f0)
{
  int index1, index2;
  for (int i = 0; i < number_of_boundaries / 2; ++i) {
    index1 = boundary_list[i * 2];
    index2 = boundary_list[i * 2 + 1];
    
    for (int j = 0; j < index1; ++j) { multi_channel_f0[i][j] = 0.0; }
    copy(f0 + index1, f0 + index2 + 1, multi_channel_f0[i] + index1);
    for (int j = index2 + 1; j < f0_length; ++j) { multi_channel_f0[i][j] = 0.0; }
  }
}

//-----------------------------------------------------------------------------
// Step 3: Voiced sections are extended based on the continuity of F0 contour
//-----------------------------------------------------------------------------
void Harvest::fixStep3(const double *f0_step2, const int f0_length,
		       const int number_of_candidates, const double * const *f0_candidates,
		       const double allowed_range, const double * const *f0_scores, double *f0_step3)
{
  copy(f0_step2, f0_step2 + f0_length, f0_step3);
  
  int *boundary_list = new int[f0_length];
  int number_of_boundaries = getBoundaryList(f0_step2, f0_length, boundary_list);

  double **multi_channel_f0 = new double *[number_of_boundaries / 2];
  for (int i = 0; i < number_of_boundaries / 2; ++i)
    { multi_channel_f0[i] = new double[f0_length]; }
  
  getMultiChannelF0(f0_step2, f0_length, boundary_list, number_of_boundaries,
		    multi_channel_f0);

  int number_of_channels =
    extend(multi_channel_f0, number_of_boundaries / 2, f0_length,
	   boundary_list, f0_candidates, number_of_candidates, allowed_range,
	   multi_channel_f0, boundary_list);
  
  mergeF0(multi_channel_f0, boundary_list, number_of_channels, f0_length,
	  f0_candidates, f0_scores, number_of_candidates, f0_step3);

  for (int i = 0; i < number_of_boundaries / 2; ++i)
    { delete[] multi_channel_f0[i]; }
  delete[] multi_channel_f0;
  delete[] boundary_list;
}

//-----------------------------------------------------------------------------
// Step 4: F0s in short unvoiced section are faked
//-----------------------------------------------------------------------------
void Harvest::fixStep4(const double *f0_step3, const int f0_length, const int threshold,
		       double *f0_step4)
{
  copy(f0_step3, f0_step3 + f0_length, f0_step4);
  
  int *boundary_list = new int[f0_length];
  int number_of_boundaries = getBoundaryList(f0_step3, f0_length, boundary_list);

  int distance;
  double tmp0, tmp1, coefficient;
  int count;
  
  for (int i = 0; i < number_of_boundaries / 2 - 1; ++i) {
    distance = boundary_list[(i + 1) * 2] - boundary_list[i * 2 + 1] - 1;
    if (distance >= threshold) continue;
    tmp0 = f0_step3[boundary_list[i * 2 + 1]] + 1;
    tmp1 = f0_step3[boundary_list[(i + 1) * 2]] - 1;
    coefficient = (tmp1 - tmp0) / (distance + 1.0);
    count = 1;
    for (int j = boundary_list[i * 2 + 1] + 1;
	 j <= boundary_list[(i + 1) * 2] - 1; ++j)
      { f0_step4[j] = tmp0 + coefficient * count++; }
  }
  
  delete[] boundary_list;
}

//-----------------------------------------------------------------------------
// FixF0Contour() obtains the likely F0 contour.
//-----------------------------------------------------------------------------
void Harvest::fixF0Contour(const double * const *f0_candidates,
			   const double * const *f0_scores,
			   const int f0_length, const int number_of_candidates,
			   double *best_f0_contour)
{
  double *tmp_f0_contour1 = new double[f0_length];
  double *tmp_f0_contour2 = new double[f0_length];

  // These parameters are optimized by speech databases.
  searchF0Base(f0_candidates, f0_scores, f0_length,
	       number_of_candidates, tmp_f0_contour1);
  fixStep1(tmp_f0_contour1, f0_length, 0.008, tmp_f0_contour2);
  fixStep2(tmp_f0_contour2, f0_length, 6, tmp_f0_contour1);
  fixStep3(tmp_f0_contour1, f0_length, number_of_candidates, f0_candidates,
	   0.18, f0_scores, tmp_f0_contour2);
  fixStep4(tmp_f0_contour2, f0_length, 9, best_f0_contour);

  delete[] tmp_f0_contour1;
  delete[] tmp_f0_contour2;
}

//-----------------------------------------------------------------------------
// This function uses zero-lag Butterworth filter.
//-----------------------------------------------------------------------------
void Harvest::filteringF0(const double *a, const double *b, double *x,
			  const int x_length, const int st, const int ed, double *y)
{
  double w[2] = { 0.0, 0.0 };
  double wt;
  double *tmp_x = new double[x_length];

  for (int i = 0; i < st; ++i) { x[i] = x[st]; }
  for (int i = ed + 1; i < x_length; ++i) { x[i] = x[ed]; }

  for (int i = 0; i < x_length; ++i) {
    wt = x[i] + a[0] * w[0] + a[1] * w[1];
    tmp_x[x_length - i - 1] = b[0] * wt + b[1] * w[0] + b[0] * w[1];
    w[1] = w[0];
    w[0] = wt;
  }

  w[0] = w[1] = 0.0;
  for (int i = 0; i < x_length; ++i) {
    wt = tmp_x[i] + a[0] * w[0] + a[1] * w[1];
    y[x_length - i - 1] = b[0] * wt + b[1] * w[0] + b[0] * w[1];
    w[1] = w[0];
    w[0] = wt;
  }

  delete[] tmp_x;
}

//-----------------------------------------------------------------------------
// SmoothF0Contour() uses the zero-lag Butterworth filter for smoothing.
//-----------------------------------------------------------------------------
void Harvest::smoothF0Contour(const double *f0, const int f0_length, double *smoothed_f0)
{
  const double b[2] =
    { 0.0078202080334971724, 0.015640416066994345 };
  const double a[2] =
    { 1.7347257688092754, -0.76600660094326412 };
  int lag = 300;
  int new_f0_length = f0_length + lag * 2;
  
  double *f0_contour = new double[new_f0_length]();
  copy(f0, f0 + f0_length, f0_contour + lag);
  
  int *boundary_list = new int[new_f0_length];
  int number_of_boundaries = getBoundaryList(f0_contour, new_f0_length, boundary_list);
  double **multi_channel_f0 = new double *[number_of_boundaries / 2];
  for (int i = 0; i < number_of_boundaries / 2; ++i)
    { multi_channel_f0[i] = new double[new_f0_length]; }
  
  getMultiChannelF0(f0_contour, new_f0_length, boundary_list,
		    number_of_boundaries, multi_channel_f0);

  for (int i = 0; i < number_of_boundaries / 2; ++i) {
    filteringF0(a, b, multi_channel_f0[i], new_f0_length,
		boundary_list[i * 2], boundary_list[i * 2 + 1], f0_contour);
    for (int j = boundary_list[i * 2]; j <= boundary_list[i * 2 + 1]; ++j)
      smoothed_f0[j - lag] = f0_contour[j];
  }

  for (int i = 0; i < number_of_boundaries / 2; ++i)
    delete[] multi_channel_f0[i];
  delete[] multi_channel_f0;
  delete[] f0_contour;
  delete[] boundary_list;
}

// static void RemoveUnreliableCandidatesSub(int i, int j,
// 					  const double * const *tmp_f0_candidates, int number_of_candidates,
// 					  double **f0_candidates, double **f0_scores) {
//   double reference_f0 = f0_candidates[i][j];
//   double error1, error2, min_error;
//   double threshold = 0.05;
//   if (reference_f0 == 0) return;
//   SelectBestF0(reference_f0, tmp_f0_candidates[i + 1],
// 	       number_of_candidates, 1.0, &error1);
//   SelectBestF0(reference_f0, tmp_f0_candidates[i - 1],
// 	       number_of_candidates, 1.0, &error2);
//   min_error = MyMinDouble(error1, error2);
//   if (min_error <= threshold) return;
//   f0_candidates[i][j] = 0;
//   f0_scores[i][j] = 0;
// }

//-----------------------------------------------------------------------------
// RemoveUnreliableCandidates().
//-----------------------------------------------------------------------------
void Harvest::removeUnreliableCandidates(const int f0_length, const int number_of_candidates,
					 double **f0_candidates, double **f0_scores)
{
  double **tmp_f0_candidates = new double *[f0_length];
  for (int i = 0; i < f0_length; ++i)
    { tmp_f0_candidates[i] = new double[number_of_candidates]; }
  
  for (int i = 1; i < f0_length - 1; ++i)
    { copy(f0_candidates[i], f0_candidates[i] + number_of_candidates, tmp_f0_candidates[i]); }
  
  for (int i = 1; i < f0_length - 1; ++i) {
    double *pnt_f0_candidates = f0_candidates[i];
    double *pnt_f0_scores = f0_scores[i];
    
    for (int j = 0; j < number_of_candidates; ++j) {
      double reference_f0 = pnt_f0_candidates[j];
      double error1, error2, min_error;
      double threshold = 0.05;
      
      if (reference_f0 == 0) continue;
      
      selectBestF0(reference_f0, tmp_f0_candidates[i + 1],
		   number_of_candidates, 1.0, error1);
      selectBestF0(reference_f0, tmp_f0_candidates[i - 1],
		   number_of_candidates, 1.0, error2);
      
      min_error = MyMinDouble(error1, error2);
     
      if (min_error <= threshold) continue;
      
      pnt_f0_candidates[j] = 0;
      pnt_f0_scores[j] = 0;
    }
  }
  
  // RemoveUnreliableCandidatesSub(i, j, tmp_f0_candidates,
  // 				number_of_candidates, f0_candidates, f0_scores);

  for (int i = 0; i < f0_length; ++i) delete[] tmp_f0_candidates[i];
  delete[] tmp_f0_candidates;
}


//-----------------------------------------------------------------------------
// GetBaseIndex() calculates the temporal positions for windowing.
//-----------------------------------------------------------------------------
void Harvest::getBaseIndex(const double current_position, const double *base_time,
			   const int base_time_length, const double fs, int *base_index)
{
  // First-aid treatment
  int basic_index = matlab_round((current_position + base_time[0]) * fs + 0.001);

  for (int i = 0; i < base_time_length; ++i) { base_index[i] = basic_index + i; }
}

//-----------------------------------------------------------------------------
// GetMainWindow() generates the window function.
//-----------------------------------------------------------------------------
void Harvest::getMainWindow(const double current_position, const int *base_index,
			    const int base_time_length, const double fs,
			    const double window_length_in_time, double *main_window)
{
  double tmp = 0.0;
  for (int i = 0; i < base_time_length; ++i) {
    tmp = (base_index[i] - 1.0) / fs - current_position;
    main_window[i] = 0.42 +
      0.5 * cos(2.0 * world::kPi * tmp / window_length_in_time) +
      0.08 * cos(4.0 * world::kPi * tmp / window_length_in_time);
  }
}

//-----------------------------------------------------------------------------
// GetDiffWindow() generates the differentiated window.
// Diff means differential.
//-----------------------------------------------------------------------------
void Harvest::getDiffWindow(const double *main_window, const int base_time_length,
			    double *diff_window)
{
  diff_window[0] = - main_window[1] / 2.0;
  diff_window[base_time_length - 1] = main_window[base_time_length - 2] / 2.0;
  
  for (int i = 1; i < base_time_length - 1; ++i)
    { diff_window[i] = -(main_window[i + 1] - main_window[i - 1]) / 2.0; }
}

//-----------------------------------------------------------------------------
// GetSpectra() calculates two spectra of the waveform windowed by windows
// (main window and diff window).
//-----------------------------------------------------------------------------
void Harvest::getSpectra(const double *x, const int x_length, const int fft_size,
			 const int *base_index, const double *main_window,
			 const double *diff_window, const int base_time_length,
			 const ForwardRealFFT *forward_real_fft, fft_complex *main_spectrum,
			 fft_complex *diff_spectrum)
{
  int *safe_index = new int[base_time_length];

  for (int i = 0; i < base_time_length; ++i)
    safe_index[i] = MyMaxInt(0, MyMinInt(x_length - 1, base_index[i] - 1));
  for (int i = 0; i < base_time_length; ++i)
    forward_real_fft->waveform[i] = x[safe_index[i]] * main_window[i];
  for (int i = base_time_length; i < fft_size; ++i)
    forward_real_fft->waveform[i] = 0.0;

  fft_execute(forward_real_fft->forward_fft);
  for (int i = 0; i <= fft_size / 2; ++i) {
    main_spectrum[i][0] = forward_real_fft->spectrum[i][0];
    main_spectrum[i][1] = -forward_real_fft->spectrum[i][1];
  }

  for (int i = 0; i < base_time_length; ++i)
    forward_real_fft->waveform[i] = x[safe_index[i]] * diff_window[i];
  for (int i = base_time_length; i < fft_size; ++i)
    forward_real_fft->waveform[i] = 0.0;
  fft_execute(forward_real_fft->forward_fft);
  for (int i = 0; i <= fft_size / 2; ++i) {
    diff_spectrum[i][0] = forward_real_fft->spectrum[i][0];
    diff_spectrum[i][1] = -forward_real_fft->spectrum[i][1];
  }

  delete[] safe_index;
}

void Harvest::fixF0(const double *power_spectrum, const double *numerator_i,
		    const int fft_size, const double fs, const double current_f0,
		    const int number_of_harmonics,
		    double *refined_f0, double *score)
{
  double *amplitude_list = new double[number_of_harmonics];
  double *instantaneous_frequency_list = new double[number_of_harmonics];

  int index;
  for (int i = 0; i < number_of_harmonics; ++i) {
    index = matlab_round(current_f0 * fft_size / fs * (i + 1));
    instantaneous_frequency_list[i] =
      (power_spectrum[index] == 0.0) ? 0.0
      : static_cast<double>(index) * fs / fft_size +
      numerator_i[index] / power_spectrum[index] * fs / 2.0 / world::kPi;
    
    amplitude_list[i] = sqrt(power_spectrum[index]);
  }
  
  double denominator = 0.0;
  double numerator = 0.0;
  *score = 0.0;
  
  for (int i = 0; i < number_of_harmonics; ++i) {
    numerator += amplitude_list[i] * instantaneous_frequency_list[i];
    denominator += amplitude_list[i] * (i + 1.0);
    *score += fabs((instantaneous_frequency_list[i] / (i + 1.0) - current_f0) / current_f0);
  }

  *refined_f0 = numerator / (denominator + world::kMySafeGuardMinimum);
  *score = 1.0 / (*score / number_of_harmonics + world::kMySafeGuardMinimum);

  delete[] amplitude_list;
  delete[] instantaneous_frequency_list;
}

//-----------------------------------------------------------------------------
// GetMeanF0() calculates the instantaneous frequency.
//-----------------------------------------------------------------------------
void Harvest::getMeanF0(const double *x, const int x_length, const double fs,
			const double current_position, const double current_f0, const int fft_size,
			const double window_length_in_time, const double *base_time,
			int base_time_length, double *refined_f0, double *refined_score)
{
  ForwardRealFFT forward_real_fft = { 0 };
  InitializeForwardRealFFT(fft_size, &forward_real_fft);
  fft_complex *main_spectrum = new fft_complex[fft_size];
  fft_complex *diff_spectrum = new fft_complex[fft_size];

  int *base_index = new int[base_time_length];
  double *main_window = new double[base_time_length];
  double *diff_window = new double[base_time_length];

  getBaseIndex(current_position, base_time, base_time_length, fs, base_index);
  
  getMainWindow(current_position, base_index, base_time_length, fs,
		window_length_in_time, main_window);
  
  getDiffWindow(main_window, base_time_length, diff_window);

  getSpectra(x, x_length, fft_size, base_index, main_window, diff_window,
	     base_time_length, &forward_real_fft, main_spectrum, diff_spectrum);

  double *power_spectrum = new double[fft_size / 2 + 1];
  double *numerator_i = new double[fft_size / 2 + 1];
  
  for (int j = 0; j <= fft_size / 2; ++j) {
    numerator_i[j] =
      main_spectrum[j][0] * diff_spectrum[j][1] - main_spectrum[j][1] * diff_spectrum[j][0];
    power_spectrum[j] =
      main_spectrum[j][0] * main_spectrum[j][0] + main_spectrum[j][1] * main_spectrum[j][1];
  }

  int number_of_harmonics = MyMinInt(static_cast<int>(fs / 2.0 / current_f0), 6);
  
  fixF0(power_spectrum, numerator_i, fft_size, fs, current_f0,
	number_of_harmonics, refined_f0, refined_score);

  

  delete[] diff_spectrum;
  delete[] diff_window;
  delete[] main_window;
  delete[] base_index;
  delete[] numerator_i;
  delete[] power_spectrum;
  delete[] main_spectrum;
  DestroyForwardRealFFT(&forward_real_fft);
}

//-----------------------------------------------------------------------------
// GetRefinedF0() calculates F0 and its score based on instantaneous frequency.
//-----------------------------------------------------------------------------
void Harvest::getRefinedF0(const double *x, const int x_length, const double fs,
			   const double current_position,
			   double current_f0,
			   double *refined_f0, double *refined_score)
{
  // double current_f0 = *refined_f0;
  
  if (current_f0 <= 0.0) {
    *refined_f0 = 0.0;
    *refined_score = 0.0;
    return;
  }

  int half_window_length = static_cast<int>(1.5 * fs / current_f0 + 1.0);
  double window_length_in_time = (2.0 * half_window_length + 1.0) / fs;
  double *base_time = new double[half_window_length * 2 + 1];
  
  for (int i = 0; i < half_window_length * 2 + 1; i++)
    base_time[i] = (-half_window_length + i) / fs;
  
  int fft_size =
    static_cast<int>(pow(2.0, 2.0 +
			 static_cast<int>(log(half_window_length * 2.0 + 1.0) / world::kLog2))
    );

  getMeanF0(x, x_length, fs, current_position, current_f0, fft_size,
	    window_length_in_time, base_time, half_window_length * 2 + 1,
	    refined_f0, refined_score);

  if (*refined_f0 < option_.f0_floor ||
      *refined_f0 > option_.f0_ceil ||
      *refined_score < 2.5) {
    *refined_f0 = 0.0;
    *refined_score = 0.0;
  }
  
  delete[] base_time;
}

//-----------------------------------------------------------------------------
// RefineF0() modifies the F0 by instantaneous frequency.
//-----------------------------------------------------------------------------
void Harvest::refineF0Candidates(const double *x, int x_length, double fs,
				 const double *temporal_positions, int f0_length, int max_candidates,
				 double **refined_f0_candidates, double **f0_scores)
{
  for (int i = 0; i < f0_length; i++) {
    for (int j = 0; j < max_candidates; ++j) {
      getRefinedF0(x, x_length, fs, temporal_positions[i],
		   refined_f0_candidates[i][j],
		   &refined_f0_candidates[i][j], &f0_scores[i][j]);
    }
  }
}

//-----------------------------------------------------------------------------
// OverlapF0Candidates() spreads the candidates to anteroposterior frames.
//-----------------------------------------------------------------------------
void Harvest::overlapF0Candidates(const int f0_length, const int number_of_candidates,
				  double **f0_candidates)
{
  int n = 3;
  for (int i = 1; i <= n; ++i)
    for (int j = 0; j < number_of_candidates; ++j) {
      
      for (int k = i; k < f0_length; ++k)
	{ f0_candidates[k][j + (number_of_candidates * i)] = f0_candidates[k - i][j]; }
      
      for (int k = 0; k < f0_length - i; ++k)
	{ f0_candidates[k][j + (number_of_candidates * (i + n))] = f0_candidates[k + i][j]; }
    }
}

//-----------------------------------------------------------------------------
// DetectOfficialF0CandidatesSub2() calculates F0 candidates in a frame
//-----------------------------------------------------------------------------
int Harvest::detectOfficialF0CandidatesSub2(const int *vuv, const double * const *raw_f0_candidates,
					    const int index, const int number_of_voiced_sections,
					    const int *st, const int *ed,
					    int max_candidates, double *f0_list)
{
  int number_of_candidates = 0;
  double tmp_f0;
  
  for (int i = 0; i < number_of_voiced_sections; ++i) {
    if (ed[i] - st[i] < 10) continue;

    tmp_f0 = 0.0;
    for (int j = st[i]; j < ed[i]; ++j)
      { tmp_f0 += raw_f0_candidates[j][index]; }
    tmp_f0 /= (ed[i] - st[i]);
    
    f0_list[number_of_candidates++] = tmp_f0;
  }

  for (int i = number_of_candidates; i < max_candidates; ++i) { f0_list[i] = 0.0; }
  
  return number_of_candidates;
}

//-----------------------------------------------------------------------------
// DetectF0CandidatesSub1() calculates VUV areas.
//-----------------------------------------------------------------------------
int Harvest::detectOfficialF0CandidatesSub1(const int *vuv,
					    const int number_of_channels,
					    int *st, int *ed)
{
  int number_of_voiced_sections = 0;
  int tmp;
  
  for (int i = 1; i < number_of_channels; ++i) {
    tmp = vuv[i] - vuv[i - 1];
    if (tmp == 1) st[number_of_voiced_sections] = i;
    if (tmp == -1) ed[number_of_voiced_sections++] = i;
  }

  return number_of_voiced_sections;
}

//-----------------------------------------------------------------------------
// DetectOfficialF0Candidates() detectes F0 candidates from multi-channel
// candidates.
//-----------------------------------------------------------------------------
int Harvest::detectOfficialF0Candidates(const double * const * raw_f0_candidates,
					const int number_of_channels, const int f0_length,
					const int max_candidates, double **f0_candidates)
{
  int number_of_candidates = 0;

  int *vuv = new int[number_of_channels];
  int *st = new int[number_of_channels];
  int *ed = new int[number_of_channels];
  int number_of_voiced_sections;
  
  for (int i = 0; i < f0_length; ++i) {
    for (int j = 0; j < number_of_channels; ++j)
      { vuv[j] = raw_f0_candidates[j][i] > 0 ? 1 : 0; }
    
    vuv[0] = vuv[number_of_channels - 1] = 0;
    
    number_of_voiced_sections =
      detectOfficialF0CandidatesSub1(vuv, number_of_channels, st, ed);
    
    number_of_candidates =
      MyMaxInt(number_of_candidates,
	       detectOfficialF0CandidatesSub2(vuv, raw_f0_candidates, i, number_of_voiced_sections,
					      st, ed, max_candidates, f0_candidates[i])
      );
  }

  delete[] vuv;
  delete[] st;
  delete[] ed;
  return number_of_candidates;
}


//-----------------------------------------------------------------------------
// CheckEvent() returns 1, provided that the input value is over 1.
// This function is for RawEventByDio().
//-----------------------------------------------------------------------------
inline int Harvest::checkEvent(int x) {
  return x > 0 ? 1 : 0;
}

//-----------------------------------------------------------------------------
// GetF0CandidateContour() calculates the F0 candidate contour in 1-ch signal.
// Calculation of F0 candidates is carried out in GetF0CandidatesSub().
//-----------------------------------------------------------------------------
void Harvest::getF0CandidateContour(const ZeroCrossings &zero_crossings, const double boundary_f0,
				    const double *temporal_positions,
				    const int f0_length, double *f0_candidate)
{
  if (0 == checkEvent(zero_crossings.number_of_negatives - 2) *
      checkEvent(zero_crossings.number_of_positives - 2) *
      checkEvent(zero_crossings.number_of_peaks - 2) *
      checkEvent(zero_crossings.number_of_dips - 2)) {
    for (int i = 0; i < f0_length; ++i) f0_candidate[i] = 0.0;
    return;
  }

  double *interpolated_f0_set[4];
  for (int i = 0; i < 4; ++i)
    interpolated_f0_set[i] = new double[f0_length];

  interp1(zero_crossings.negative_interval_locations,
	  zero_crossings.negative_intervals,
	  zero_crossings.number_of_negatives,
	  temporal_positions, f0_length, interpolated_f0_set[0]);
  interp1(zero_crossings.positive_interval_locations,
	  zero_crossings.positive_intervals,
	  zero_crossings.number_of_positives,
	  temporal_positions, f0_length, interpolated_f0_set[1]);
  interp1(zero_crossings.peak_interval_locations,
	  zero_crossings.peak_intervals, zero_crossings.number_of_peaks,
	  temporal_positions, f0_length, interpolated_f0_set[2]);
  interp1(zero_crossings.dip_interval_locations,
	  zero_crossings.dip_intervals, zero_crossings.number_of_dips,
	  temporal_positions, f0_length, interpolated_f0_set[3]);

  // GetF0CandidateContourSub(interpolated_f0_set, f0_length, boundary_f0, f0_candidate);
  // GetF0CandidateContourSub
  double upper = boundary_f0 * 1.1;
  double lower = boundary_f0 * 0.9;
  
  for (int i = 0; i < f0_length; ++i) {
    f0_candidate[i] = (interpolated_f0_set[0][i] +
		       interpolated_f0_set[1][i] + interpolated_f0_set[2][i] +
		       interpolated_f0_set[3][i]) / 4.0;

    if (f0_candidate[i] > upper || f0_candidate[i] < lower ||
	f0_candidate[i] > option_.f0_ceil || f0_candidate[i] < option_.f0_floor)
      f0_candidate[i] = 0.0;
  }
  
  for (int i = 0; i < 4; ++i) delete[] interpolated_f0_set[i];
}

//-----------------------------------------------------------------------------
// struct for RawEventByHarvest()
// "negative" means "zero-crossing point going from positive to negative"
// "positive" means "zero-crossing point going from negative to positive"
//-----------------------------------------------------------------------------

void Harvest::ZeroCrossings::init_array(const int n)
{
  negative_interval_locations = new double[n];
  negative_intervals = new double[n];
  positive_interval_locations = new double[n];
  positive_intervals = new double[n];
  peak_interval_locations = new double[n];
  peak_intervals = new double[n];
  dip_interval_locations = new double[n];
  dip_intervals = new double[n];
}

void Harvest::ZeroCrossings::destroy()
{
  delete[] negative_interval_locations;
  delete[] positive_interval_locations;
  delete[] peak_interval_locations;
  delete[] dip_interval_locations;
  delete[] negative_intervals;
  delete[] positive_intervals;
  delete[] peak_intervals;
  delete[] dip_intervals;
}

// typedef static struct {
//   double *negative_interval_locations;
//   double *negative_intervals;
//   int number_of_negatives;
//   double *positive_interval_locations;
//   double *positive_intervals;
//   int number_of_positives;
//   double *peak_interval_locations;
//   double *peak_intervals;
//   int number_of_peaks;
//   double *dip_interval_locations;
//   double *dip_intervals;
//   int number_of_dips;

//   void init_array(const int n)
//   {
//     negative_interval_locations = new double[n];
//     negative_intervals = new double[n];
//     positive_interval_locations = new double[n];
//     positive_intervals = new double[n];
//     peak_interval_locations = new double[n];
//     peak_intervals = new double[n];
//     dip_interval_locations = new double[n];
//     dip_intervals = new double[n];
//   }
  
//   void destroy()
//   {
//     delete[] negative_interval_locations;
//     delete[] positive_interval_locations;
//     delete[] peak_interval_locations;
//     delete[] dip_interval_locations;
//     delete[] negative_intervals;
//     delete[] positive_intervals;
//     delete[] peak_intervals;
//     delete[] dip_intervals;
//   }
// } ZeroCrossings;

//-----------------------------------------------------------------------------
// ZeroCrossingEngine() calculates the zero crossing points from positive to
// negative.
//-----------------------------------------------------------------------------
int Harvest::zeroCrossingEngine(const double *filtered_signal, const int y_length,
				const double fs, double *interval_locations, double *intervals)
{
  int *negative_going_points = new int[y_length];

  for (int i = 0; i < y_length - 1; ++i)
    { negative_going_points[i] =
	(0.0 < filtered_signal[i] && filtered_signal[i + 1] <= 0.0) ? i + 1 : 0; }
  
  negative_going_points[y_length - 1] = 0;

  int *edges = new int[y_length];
  int count = 0;
  
  for (int i = 0; i < y_length; ++i)
    { if (negative_going_points[i] > 0)
	{ edges[count++] = negative_going_points[i]; } }

  if (count < 2) {
    delete[] edges;
    delete[] negative_going_points;
    return 0;
  }

  double *fine_edges = new double[count];
  
  for (int i = 0; i < count; ++i)
    { fine_edges[i] = edges[i] - filtered_signal[edges[i] - 1] /
	(filtered_signal[edges[i]] - filtered_signal[edges[i] - 1]); }

  for (int i = 0; i < count - 1; ++i) {
    intervals[i] = fs / (fine_edges[i + 1] - fine_edges[i]);
    interval_locations[i] = (fine_edges[i] + fine_edges[i + 1]) / 2.0 / fs;
  }

  delete[] fine_edges;
  delete[] edges;
  delete[] negative_going_points;
  return count - 1;
}

//-----------------------------------------------------------------------------
// GetFourZeroCrossingIntervals() calculates four zero-crossing intervals.
// (1) Zero-crossing going from negative to positive.
// (2) Zero-crossing going from positive to negative.
// (3) Peak, and (4) dip. (3) and (4) are calculated from the zero-crossings of
// the differential of waveform.
//-----------------------------------------------------------------------------
void Harvest::getFourZeroCrossingIntervals(double *filtered_signal, const int y_length,
					   const double actual_fs, ZeroCrossings &zero_crossings)
{
  int maximum_number = y_length;
  zero_crossings.init_array(maximum_number);

  zero_crossings.number_of_negatives =
    zeroCrossingEngine(filtered_signal, y_length, actual_fs,
		       zero_crossings.negative_interval_locations,
		       zero_crossings.negative_intervals);

  for (int i = 0; i < y_length; ++i) { filtered_signal[i] *= -1; }
  
  zero_crossings.number_of_positives =
    zeroCrossingEngine(filtered_signal, y_length, actual_fs,
		       zero_crossings.positive_interval_locations,
		       zero_crossings.positive_intervals);

  for (int i = 0; i < y_length - 1; ++i) { filtered_signal[i] -= filtered_signal[i + 1]; }
  
  zero_crossings.number_of_peaks =
    zeroCrossingEngine(filtered_signal, y_length - 1, actual_fs,
		       zero_crossings.peak_interval_locations,
		       zero_crossings.peak_intervals);
  
  for (int i = 0; i < y_length - 1; ++i) { filtered_signal[i] *= -1; }
  
  zero_crossings.number_of_dips =
    zeroCrossingEngine(filtered_signal, y_length - 1, actual_fs,
		       zero_crossings.dip_interval_locations,
		       zero_crossings.dip_intervals);
}

//-----------------------------------------------------------------------------
// GetFilteredSignal() calculates the signal that is the convolution of the
// input signal and band-pass filter.
//-----------------------------------------------------------------------------
void Harvest::getFilteredSignal(const double boundary_f0, const int fft_size, const double fs,
				const fft_complex *y_spectrum, const int y_length,
				double *filtered_signal)
{  
  int filter_length_half = matlab_round(fs / boundary_f0 * 2.0);
  double *band_pass_filter = new double[fft_size]();
  NuttallWindow(filter_length_half * 2 + 1, band_pass_filter);
  
  for (int i = -filter_length_half; i <= filter_length_half; ++i)
    { band_pass_filter[i + filter_length_half] *= cos(2 * world::kPi * boundary_f0 * i / fs); }
  
  // fft_complex *band_pass_filter_spectrum = new fft_complex[fft_size];
  fft_complex *band_pass_filter_spectrum = new fft_complex[fft_size / 2 + 1];
  fft_plan forwardFFT = fft_plan_dft_r2c_1d(fft_size, band_pass_filter,
					    band_pass_filter_spectrum, FFT_ESTIMATE);
  fft_execute(forwardFFT);

  // Convolution
  double tmp = y_spectrum[0][0] * band_pass_filter_spectrum[0][0] -
    y_spectrum[0][1] * band_pass_filter_spectrum[0][1];
  band_pass_filter_spectrum[0][1] =
    y_spectrum[0][0] * band_pass_filter_spectrum[0][1] +
    y_spectrum[0][1] * band_pass_filter_spectrum[0][0];
  band_pass_filter_spectrum[0][0] = tmp;
  
  for (int i = 1; i <= fft_size / 2; ++i) {
    tmp = y_spectrum[i][0] * band_pass_filter_spectrum[i][0] -
      y_spectrum[i][1] * band_pass_filter_spectrum[i][1];
    band_pass_filter_spectrum[i][1] =
      y_spectrum[i][0] * band_pass_filter_spectrum[i][1] +
      y_spectrum[i][1] * band_pass_filter_spectrum[i][0];
    band_pass_filter_spectrum[i][0] = tmp;
    // band_pass_filter_spectrum[fft_size - i - 1][0] =
    //   band_pass_filter_spectrum[i][0];
    // band_pass_filter_spectrum[fft_size - i - 1][1] =
    //   band_pass_filter_spectrum[i][1];
  }

  fft_plan inverseFFT = fft_plan_dft_c2r_1d(fft_size,
					    band_pass_filter_spectrum, filtered_signal, FFT_ESTIMATE);
  fft_execute(inverseFFT);

  // Compensation of the delay.
  int index_bias = filter_length_half + 1;
  rotate(filtered_signal, filtered_signal + index_bias, filtered_signal + fft_size);
  // for (int i = 0; i < y_length; ++i)
  //   filtered_signal[i] = filtered_signal[i + index_bias];

  fft_destroy_plan(inverseFFT);
  fft_destroy_plan(forwardFFT);
  delete[] band_pass_filter_spectrum;
  delete[] band_pass_filter;
}

//-----------------------------------------------------------------------------
// GetF0CandidateFromRawEvent() f0 candidate contour in 1-ch signal
//-----------------------------------------------------------------------------
void Harvest::getF0CandidateFromRawEvent(const double boundary_f0, const double fs,
					 const fft_complex *y_spectrum, const int y_length, const int fft_size,
					 const double *temporal_positions, const int f0_length,
					 double *f0_candidate) {
  double *filtered_signal = new double[fft_size];
  getFilteredSignal(boundary_f0, fft_size, fs, y_spectrum,
		    y_length, filtered_signal);

  ZeroCrossings zero_crossings;
  getFourZeroCrossingIntervals(filtered_signal, y_length, fs, zero_crossings);

  getF0CandidateContour(zero_crossings, boundary_f0, temporal_positions, f0_length, f0_candidate);

  // DestroyZeroCrossings(&zero_crossings);
  zero_crossings.destroy();
  delete[] filtered_signal;
}

//-----------------------------------------------------------------------------
// GetRawF0Candidates() calculates f0 candidates in all channels.
//-----------------------------------------------------------------------------
void Harvest::getRawF0Candidates(const double *boundary_f0_list,
				 const int number_of_bands, const double actual_fs, const int y_length,
				 const double *temporal_positions, const int f0_length,
				 const fft_complex *y_spectrum, const int fft_size,
				 double **raw_f0_candidates) {
  for (int i = 0; i < number_of_bands; ++i) {
    { getF0CandidateFromRawEvent(boundary_f0_list[i], actual_fs, y_spectrum,
				 y_length, fft_size, temporal_positions, f0_length,
				 raw_f0_candidates[i]); }
  }
}

//-----------------------------------------------------------------------------
// HarvestGeneralBodySub() is the subfunction of HarvestGeneralBody()
//-----------------------------------------------------------------------------
int Harvest::generalBodySub(const double *boundary_f0_list,
			    const int number_of_channels, const int f0_length,
			    const double actual_fs, const int y_length,
			    const double *temporal_positions, const fft_complex *y_spectrum,
			    const int fft_size, const int max_candidates, double **f0_candidates)
{
  double **raw_f0_candidates = new double *[number_of_channels];
  for (int i = 0; i < number_of_channels; ++i)
    raw_f0_candidates[i] = new double[f0_length];
  
  getRawF0Candidates(boundary_f0_list, number_of_channels,
		     actual_fs, y_length, temporal_positions, f0_length, y_spectrum,
		     fft_size, raw_f0_candidates);
  
  int number_of_candidates =
    detectOfficialF0Candidates(raw_f0_candidates,
			       number_of_channels, f0_length, max_candidates, f0_candidates);

  overlapF0Candidates(f0_length, number_of_candidates, f0_candidates);
  
  for (int i = 0; i < number_of_channels; ++i)
    delete[] raw_f0_candidates[i];
  delete[] raw_f0_candidates;
  
  return number_of_candidates;
}


//-----------------------------------------------------------------------------
// HarvestGeneralBody() estimates the F0 contour based on Harvest.
//-----------------------------------------------------------------------------
void Harvest::generalBody(const double *x, const int x_length, const int fs,
			  const int frame_period, const double channels_in_octave,
			  const int speed,
			  double *temporal_positions, double *f0)
{
  double adjusted_f0_floor = option_.f0_floor * 0.9;
  double adjusted_f0_ceil = option_.f0_ceil * 1.1;
  int number_of_channels = 1 + static_cast<int>(log(adjusted_f0_ceil / adjusted_f0_floor) /
						world::kLog2 * channels_in_octave);
  double *boundary_f0_list = new double[number_of_channels];
  
  for (int i = 0; i < number_of_channels; ++i)
    { boundary_f0_list[i] = adjusted_f0_floor * pow(2.0, (i + 1) / channels_in_octave); }

  // normalization
  int decimation_ratio = MyMaxInt(MyMinInt(speed, 12), 1);
  int y_length = (1 + static_cast<int>(x_length / decimation_ratio));
  double actual_fs = static_cast<double>(fs) / decimation_ratio;
  int fft_size =
    GetSuitableFFTSize(y_length +
		       (4 * static_cast<int>(1.0 + actual_fs / boundary_f0_list[0] / 2.0)));
  
  // Calculation of the spectrum used for the f0 estimation
  double *y = new double[fft_size](); // init
  fft_complex *y_spectrum = new fft_complex[fft_size / 2 + 1];
  
  getWaveformAndSpectrum(x, x_length, y_length, actual_fs, fft_size,
			 decimation_ratio, y, y_spectrum);

  int f0_length = getSamples(fs, x_length, frame_period);
  
  for (int i = 0; i < f0_length; ++i) {
    temporal_positions[i] = i * frame_period / 1000.0;
    f0[i] = 0.0;
  }

  int overlap_parameter = 7;
  int max_candidates =
    matlab_round(number_of_channels / 10) * overlap_parameter;
  double **f0_candidates = new double *[f0_length];
  double **f0_candidates_score = new double *[f0_length];
  
  for (int i = 0; i < f0_length; ++i) {
    f0_candidates[i] = new double[max_candidates]();
    f0_candidates_score[i] = new double[max_candidates]();
  }

  DWORD t4 = timeGetTime();

  int number_of_candidates =
    generalBodySub(boundary_f0_list,
		   number_of_channels, f0_length, actual_fs, y_length, temporal_positions,
		   y_spectrum, fft_size, max_candidates, f0_candidates)
    * overlap_parameter;
  
  DWORD t5 = timeGetTime();
  
  refineF0Candidates(y, y_length, actual_fs, temporal_positions, f0_length,
		     number_of_candidates, f0_candidates, f0_candidates_score);
  
  removeUnreliableCandidates(f0_length, number_of_candidates,
			     f0_candidates, f0_candidates_score);

  double *best_f0_contour = new double[f0_length];
  
  fixF0Contour(f0_candidates, f0_candidates_score, f0_length,
	       number_of_candidates, best_f0_contour);
  
  smoothF0Contour(best_f0_contour, f0_length, f0);

  DWORD t6 = timeGetTime();

  delete[] y;
  delete[] best_f0_contour;
  delete[] y_spectrum;
  for (int i = 0; i < f0_length; ++i) {
    delete[] f0_candidates[i];
    delete[] f0_candidates_score[i];
  }
  delete[] f0_candidates;
  delete[] f0_candidates_score;
  delete[] boundary_f0_list;
  
  printf("Harvest time4: %d [msec]\n", t5 - t4);
  printf("Harvest time5: %d [msec]\n", t6 - t5);
}

} // end namespace world_class
