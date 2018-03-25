//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise,
// Copyright 2018 Yukara Ikemiya
//-----------------------------------------------------------------------------

#ifndef WORLD_CLASS_HARVEST_HPP
#define WORLD_CLASS_HARVEST_HPP

#include "world_fft.hpp"
#include "world_common.hpp"


namespace world_class
{

typedef struct HarvestOption{
	double f0_floor;
	double f0_ceil;
	double frame_period;

	double target_fs;
	double channels_in_octave;

	bool use_cos_table;
  
	HarvestOption();
	void copy(const HarvestOption& option);
} HarvestOption;


class Harvest
{

public:
	
	Harvest(const int fs, const HarvestOption &option);
	~Harvest();

	void compute(
		const double* x, int x_length,double *temporal_positions, double *f0
	);

	int getSamples(int fs, int x_length, double frame_period);
	int getSamples(int fs, int x_length);

private:
  
	HarvestOption option_;
	int fs_;
	int decimation_ratio_;
	double actual_fs_;
	int num_thread_;
  
	const double *x_;
	int x_length_;
	double *y_;
	int y_length_;
	double *temporal_positions_;
	
	int f0_length_;
	int number_of_candidates_;
  
	double **f0_candidates_;
	double **f0_candidates_score_;

	// buffers for getMeanF0
	int max_fft_size_;
	int max_base_time_length_;
	fft_complex *main_spectrum_;
	fft_complex *diff_spectrum_;
	int *base_index_;
	double *main_window_;
	double *diff_window_;
	double *power_spectrum_;
	double *numerator_i_;
	// buffers for refineF0Candidates
	double *base_time_;
	// buffers for getSpectra
	int *safe_index_;

	// FFTs
	ForwardRealFFT *structFFTs_;
	int first_fft_index_;
	int num_fft_;

	// cosine table
	int num_cos_div_;
	double *cos_table_;

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

		void initialize(int n);
		void destroy();
	};

	void generalBody(
		const double *x, int x_length,
		int frame_period, double channels_in_octave,
		double *temporal_positions, double *f0
	);
	int generalBodySub(
		const double *boundary_f0_list, int number_of_channels,
		const fft_complex *y_spectrum, int fft_size, int max_candidates
	);

	void prepareFFTs();
	void destroyFFTs();

	void get_cos_table();

	void getWaveformAndSpectrum(
		int fft_size, int decimation_ratio, fft_complex *y_spectrum
	);

	void getRawF0Candidates(
		const double *boundary_f0_list, int number_of_bands,
		const fft_complex *y_spectrum, int fft_size,
		double **raw_f0_candidates
	);
  
	void getFilteredSignal(
		double boundary_f0, int fft_size, const fft_complex *y_spectrum,
		double *filtered_signal
	);

	void getFourZeroCrossingIntervals(
		double *filtered_signal, ZeroCrossings &zero_crossings
	);
	int zeroCrossingEngine(
		const double *filtered_signal, int y_length,
		double fs, double *interval_locations, double *intervals
	);

	void getF0CandidateContour(
		const ZeroCrossings &zero_crossings, double boundary_f0,
		double *f0_candidate
	);
	inline int checkEvent(int x);

	int detectOfficialF0Candidates(
		const double * const * raw_f0_candidates,
		int number_of_channels, int f0_length, int max_candidates,
		double **f0_candidates
	);
	int detectOfficialF0CandidatesSub1(
		const int *vuv, int number_of_channels,
		int *st, int *ed
	);
	int detectOfficialF0CandidatesSub2(
		const double * const *raw_f0_candidates,
		int index, int number_of_voiced_sections,
		const int *st, const int *ed,
		int max_candidates, double *f0_list
	);

	void overlapF0Candidates(
		int f0_length, int number_of_candidates,
		double **f0_candidates
	);

	void refineF0Candidates();
	void getMeanF0(
		double current_position, double current_f0, int fft_index,
		double window_length_in_time, const double *base_time,
		int base_time_length, double *refined_f0, double *refined_score,
		int thread_id = 0
	);
	void fixF0(
		const double *power_spectrum, const double *numerator_i,
		int fft_size, double fs, double current_f0, int number_of_harmonics,
		double *refined_f0, double *score
	);
	void getSpectra(
		const double *x, int x_length, int fft_size,
		const int *base_index, const double *main_window,
		const double *diff_window, int base_time_length,
		const ForwardRealFFT &forward_real_fft, fft_complex *main_spectrum,
		fft_complex *diff_spectrum, int thread_id = 0
	);
	inline void getDiffWindow(
		const double *main_window, int base_time_length,
		double *diff_window
	);
	inline void getMainWindow(
		double current_position, const int *base_index,
		int base_time_length, double fs, double window_length_in_time,
		double *main_window
	);
	void getBaseIndex(
		double current_position, const double *base_time,
		int base_time_length, double fs, int *base_index
	);

	void removeUnreliableCandidates();
	double selectBestF0(
		double reference_f0, const double *f0_candidates,
		int number_of_candidates, double allowed_range, double &best_error
	);

	void smoothF0Contour(const double *f0, double *smoothed_f0);
	void filteringF0(
		const double *a, const double *b, double *x,
		int x_length, int st, int ed, double *y
	);

	void fixF0Contour(double *best_f0_contour);
	void fixStep4(const double *f0_step3, int threshold, double *f0_step4);
	void fixStep3(const double *f0_step2, double allowed_range, double *f0_step3);
	void fixStep2(const double *f0_step1, int voice_range_minimum, double *f0_step2);
	void fixStep1(const double *f0_base, double allowed_range, double *f0_step1);
	void searchF0Base(
		const double * const *f0_candidates, const double * const *f0_scores,
		int f0_length, int number_of_candidates, double *base_f0_contour
	);

	int getBoundaryList(const double *f0, int f0_length, int *boundary_list);
	void mergeF0(
		const double * const *multi_channel_f0, int *boundary_list,
		int number_of_channels, int f0_length, const double * const *f0_candidates,
		const double * const *f0_scores, int number_of_candidates,
		double *merged_f0
	);
	int mergeF0Sub(
		double *merged_f0, int st1, int ed1,
		const double *f0_2, int st2, int ed2, const double * const *f0_candidates,
		const double * const *f0_scores, int number_of_candidates
	);
	inline double searchScore(
		double f0, const double *f0_candidates,
		const double  *f0_scores, int number_of_candidates
	);
	int extend(
		double **extended_f0, int number_of_sections, int f0_length,
		int *boundary_list, const double * const *f0_candidates,
		int number_of_candidates, double allowed_range
	);
	int extendF0(
		double *extended_f0, int origin,
		int last_point, int shift, const double * const *f0_candidates,
		int number_of_candidates, double allowed_range
	);
	void getMultiChannelF0(
		const double *f0, int f0_length,
		const int *boundary_list, int number_of_boundaries,
		double **multi_channel_f0
	);
	inline void swapArray(int index1, int index2, double **f0, int *boundary);
	inline int MyAbsInt(int x);

  
};

} // end namespace world_class

#endif
