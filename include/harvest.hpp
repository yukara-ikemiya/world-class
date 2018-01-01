#ifndef WORLD_CLASS_HARVEST_HPP
#define WORLD_CLASS_HARVEST_HPP

#include "world_fft.h"
#include "world_common.h"


namespace world_class
{

typedef struct HarvestOption{
  double f0_floor;
  double f0_ceil;
  double frame_period;
  HarvestOption();
} HarvestOption;


class Harvest
{

public:

  Harvest();
  Harvest(const HarvestOption& option);

  void compute(const double* x, int x_length, int fs,
	       double *temporal_positions, double *f0);

  int getSamples(int fs, int x_length, double frame_period);
  int getSamples(int fs, int x_length);

private:

  struct ZeroCrossings{
    double *negative_interval_locations;
    double *negative_intervals;
    int number_of_negatives;
    double *positive_interval_locations;
    double *positive_intervals;
    int number_of_positives;
    double *peak_interval_locations;
    double *peak_intervals;
    int number_of_peaks;
    double *dip_interval_locations;
    double *dip_intervals;
    int number_of_dips;

    void init_array(int n);
    void destroy();
  };
  
  HarvestOption option_;

  void generalBody(const double *x, int x_length, int fs,
		   int frame_period, double channels_in_octave,
		   int speed, double *temporal_positions, double *f0);
  int generalBodySub(const double *boundary_f0_list,
		     int number_of_channels, int f0_length, double actual_fs, int y_length,
		     const double *temporal_positions, const fft_complex *y_spectrum,
		     int fft_size, int max_candidates, double **f0_candidates);

  void getWaveformAndSpectrum(const double *x, int x_length, int y_length, double actual_fs,
			      int fft_size, int decimation_ratio,
			      double *y, fft_complex *y_spectrum);
  void getWaveformAndSpectrumSub(const double *x, int x_length, int y_length, double actual_fs,
				 int decimation_ratio, double *y);

  void getRawF0Candidates(const double *boundary_f0_list,
			  int number_of_bands, double actual_fs, int y_length,
			  const double *temporal_positions, int f0_length,
			  const fft_complex *y_spectrum, int fft_size,
			  double **raw_f0_candidates);
  void getF0CandidateFromRawEvent(double boundary_f0, double fs,
				  const fft_complex *y_spectrum, int y_length, int fft_size,
				  const double *temporal_positions, int f0_length,
				  double *f0_candidate);

  void getFilteredSignal(double boundary_f0, int fft_size, double fs,
			 const fft_complex *y_spectrum, int y_length,
			 double *filtered_signal);

  void getFourZeroCrossingIntervals(double *filtered_signal, int y_length,
				    double actual_fs, ZeroCrossings &zero_crossings);
  int zeroCrossingEngine(const double *filtered_signal, int y_length,
			 double fs, double *interval_locations, double *intervals);

  void getF0CandidateContour(const ZeroCrossings &zero_crossings, double boundary_f0,
			     const double *temporal_positions,
			     int f0_length, double *f0_candidate);
  inline int checkEvent(int x);

  int detectOfficialF0Candidates(const double * const * raw_f0_candidates,
				 int number_of_channels, int f0_length, int max_candidates,
				 double **f0_candidates);
  int detectOfficialF0CandidatesSub1(const int *vuv,
				     int number_of_channels,
				     int *st, int *ed);
  int detectOfficialF0CandidatesSub2(const int *vuv, const double * const *raw_f0_candidates,
				     int index, int number_of_voiced_sections,
				     const int *st, const int *ed,
				     int max_candidates, double *f0_list);

  void overlapF0Candidates(int f0_length, int number_of_candidates,
			   double **f0_candidates);

  void refineF0Candidates(const double *x, int x_length, double fs,
			  const double *temporal_positions, int f0_length, int max_candidates,
			  double **refined_f0_candidates, double **f0_scores);
  void getRefinedF0(const double *x, int x_length, double fs, double current_position,
		    double current_f0, double *refined_f0, double *refined_score);
  void getMeanF0(const double *x, int x_length, double fs,
		 double current_position, double current_f0, int fft_size,
		 double window_length_in_time, const double *base_time,
		 int base_time_length, double *refined_f0, double *refined_score);
  void fixF0(const double *power_spectrum, const double *numerator_i,
	     int fft_size, double fs, double current_f0, int number_of_harmonics,
	     double *refined_f0, double *score);
  void getSpectra(const double *x, int x_length, int fft_size,
		  const int *base_index, const double *main_window,
		  const double *diff_window, int base_time_length,
		  const ForwardRealFFT *forward_real_fft, fft_complex *main_spectrum,
		  fft_complex *diff_spectrum);
  void getDiffWindow(const double *main_window, int base_time_length,
		     double *diff_window);
  void getMainWindow(double current_position, const int *base_index,
		     int base_time_length, double fs, double window_length_in_time,
		     double *main_window);
  void getBaseIndex(double current_position, const double *base_time,
		    int base_time_length, double fs, int *base_index);

  void removeUnreliableCandidates(int f0_length, int number_of_candidates,
				  double **f0_candidates, double **f0_scores);
  double selectBestF0(double reference_f0, const double *f0_candidates,
		      int number_of_candidates, double allowed_range, double &best_error);

  void smoothF0Contour(const double *f0, int f0_length, double *smoothed_f0);
  void filteringF0(const double *a, const double *b, double *x,
		   int x_length, int st, int ed, double *y);

  void fixF0Contour(const double * const *f0_candidates, const double * const *f0_scores,
		    int f0_length, int number_of_candidates, double *best_f0_contour);
  void fixStep4(const double *f0_step3, int f0_length, int threshold, double *f0_step4);
  void fixStep3(const double *f0_step2, int f0_length,
		int number_of_candidates, const double * const *f0_candidates,
		double allowed_range, const double * const *f0_scores, double *f0_step3);
  void fixStep2(const double *f0_step1, int f0_length,
		int voice_range_minimum, double *f0_step2);
  void fixStep1(const double *f0_base, int f0_length,
		double allowed_range, double *f0_step1);
  void searchF0Base(const double * const *f0_candidates, const double * const *f0_scores,
		    int f0_length, int number_of_candidates, double *base_f0_contour);

  int getBoundaryList(const double *f0, int f0_length, int *boundary_list);
  void mergeF0(const double * const *multi_channel_f0, int *boundary_list,
	       int number_of_channels, int f0_length, const double * const *f0_candidates,
	       const double * const *f0_scores, int number_of_candidates,
	       double *merged_f0);
  int mergeF0Sub(const double *f0_1, int f0_length, int st1, int ed1,
		 const double *f0_2, int st2, int ed2, const double * const *f0_candidates,
		 const double * const *f0_scores, int number_of_candidates,
		 double *merged_f0);
  inline double searchScore(double f0, const double *f0_candidates,
			    const double  *f0_scores, int number_of_candidates);
  int extend(const double * const *multi_channel_f0,
	     int number_of_sections, int f0_length, const int *boundary_list,
	     const double * const *f0_candidates, int number_of_candidates,
	     double allowed_range, double **extended_f0, int *shifted_boundary_list);
  int extendF0(const double *f0, int f0_length, int origin,
	       int last_point, int shift, const double * const *f0_candidates,
	       int number_of_candidates, double allowed_range, double *extended_f0);
  void getMultiChannelF0(const double *f0, int f0_length,
			 const int *boundary_list, int number_of_boundaries,
			 double **multi_channel_f0);
  inline void swapArray(int index1, int index2, double **f0, int *boundary);
  inline int MyAbsInt(int x);

  
};

} // end namespace world_class

#endif
