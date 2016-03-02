#include <R.h>
#include <Rdefines.h>

// math.c
int find_max(int num_data, int* data);
void gaussian_kernel(double* kernel,int size, double standard_deviation);
void gaussian_kernel_2d(double* kernel,int x_size,int y_size,double standard_deviation);
void smooth_double_gaussian(double* array, int array_size, double smooth, double invalid);
SEXP smooth_double_gaussian_cwrap(SEXP array_r, SEXP array_size_r, SEXP smooth_r, SEXP invalid_r);
void smooth_double_gaussian_degrees(double* array, int array_size,double smooth, double invalid);
SEXP smooth_double_gaussian_degrees_cwrap(SEXP array_r, SEXP array_size_r, SEXP sd_r, SEXP invalid_r);
void set_array_to_value_int (int* array, int array_size, int value);
void set_array_to_value_double (double* array, int array_size, double value);
void smooth_double_gaussian_2d(double* array, int x_size,int y_size, double smooth, double invalid);
SEXP smooth_double_gaussian_degrees_cwrap(SEXP array_r, SEXP array_size_r, SEXP sd_r, SEXP invalid_r);


// spatial.c
SEXP speed_from_whl_cwrap(SEXP x_whl_r,
			  SEXP y_whl_r,
			  SEXP whl_lines_r,
			  SEXP look_back_max_r,
			  SEXP look_ahead_max_r,
			  SEXP px_per_cm_r,
			  SEXP res_sampling_rate_r, // ex 20 000
			  SEXP res_samples_per_whl_sample_r); // ex 512
void speed_from_whl(double* x_whl,
		    double* y_whl,
		    double* speed,
		    int whl_lines,
		    int look_back_max,
		    int look_ahead_max,
		    double px_per_cm,
		    int res_sampling_rate, // ex 20 000
		    int res_samples_per_whl_sample);
SEXP angular_speed_from_hd_cwrap(SEXP hd_r,
				 SEXP whl_lines_r,
				 SEXP look_back_max_r,
				 SEXP look_ahead_max_r,
				 SEXP res_sampling_rate_r, 
				 SEXP res_samples_per_whl_sample_r);
void angular_speed_from_hd(double* hd, 
			   double* angular_speed, 
			   int whl_lines, 
			   int look_back_max, 
			   int look_ahead_max, 
			   int res_sampling_rate, 
			   int res_samples_per_whl_sample);

SEXP speed_intervals_cwrap(SEXP speed_r, SEXP whl_lines_r, SEXP res_samples_per_whl_sample_r, SEXP min_speed_r, SEXP max_speed_r);
int speed_intervals_count(double* speed, int whl_lines, int res_samples_per_whl_sample,double min_speed,double max_speed);
void speed_intervals(double* speed, int whl_lines, int res_samples_per_whl_sample,double min_speed,double max_speed,int* start,int* end);

void spike_position_no_interval(double *x_whl,int whl_lines,int *res, int res_lines, double *x_spike, int res_samples_per_whl_sample);
SEXP speed_at_res_values_cwrap(SEXP speed_r, SEXP whl_lines_r, SEXP res_r, SEXP res_lines_r, SEXP res_samples_per_whl_sample_r);
SEXP spike_position_cwrap(SEXP x_whl_r, SEXP y_whl_r, SEXP whl_lines_r, SEXP res_r, SEXP res_lines_r, SEXP res_samples_per_whl_sample_r, SEXP start_interval_r, SEXP end_interval_r, SEXP interval_lines_r);
void spike_position(double *x_whl, double *y_whl,int whl_lines,int *res,int res_lines, double *x_spike, double *y_spike, int res_samples_per_whl_sample, int* start_interval, int* end_interval, int interval_lines);
void occupancy_map(int x_bins, int y_bins,  double pixels_per_bin_x, double pixels_per_bin_y, double *x_whl, double *y_whl, int whl_lines, double *map, double ms_per_sample, int *start_interval,  int *end_interval, int interval_lines,  int res_samples_per_whl_sample);
SEXP occupancy_map_cwrap(SEXP x_bins_r, SEXP y_bins_r, SEXP pixels_per_bin_x_r, SEXP pixels_per_bin_y_r, SEXP x_whl_r,  SEXP y_whl_r, SEXP whl_lines_r,  SEXP ms_per_sample_r, SEXP start_interval_r, SEXP end_interval_r,  SEXP interval_lines_r, SEXP res_samples_per_whl_sample_r);
void create_place_field( int x_bins, int y_bins, double pixels_per_bin_x, double pixels_per_bin_y, double *x_spike,  double *y_spike, int *clu, int res_lines, int target_cell, double *occupancy_map, double *place_field);
SEXP firing_rate_map_2d_cwrap(SEXP num_bins_x_r, SEXP num_bins_y_r, SEXP pixels_per_bin_x_r,SEXP pixels_per_bin_y_r, SEXP x_spike_r, SEXP y_spike_r, SEXP clu_r,SEXP res_lines_r, SEXP cells_r, SEXP num_cells_r, SEXP occ_map_r, SEXP smooth_map_sd_r);




// interval.c
int check_interval_chronology_between(int num_lines, 
				      int* start, 
				      int* end);

void set_res_outside_interval_to_minus_one(int interval_lines,
					   int* start_interval_index,
					   int* end_interval_index,
					   int res_lines,
					   int* res);
// crosscorrelation.c
void crosscorrelation_one_cell(int clu_no1, // cell of interest1
			       int clu_no2, // cell of interest2
			       int* clu,  
			       int* res, 
			       int res_lines,
			       int* histo, // pointer to one histogram
			       int histo_size, // histo_size
			       int window_size, // size of the window_size in res value
			       int* start_interval_index,
			       int* end_interval_index,
			       int interval_lines);
void crosscorrelation_one_cell_probability(int clu_no1, // cell of interest1
					   int clu_no2, // cell of interest2
					   int* clu,  
					   int* res, 
					   int res_lines,
					   double* histo, // pointer to one histogram
					   int histo_size, // histo_size
					   int window_size, // size of the window_size in res value
					   int* start_interval_index,
					   int* end_interval_index,
					   int interval_lines);
void crosscorrelation(int* clu1, // cells of interest
		      int* clu2,
		      int cell_pair_lines,
		      int* clu,  
		      int* res, 
		      int res_lines,
		      int* histo, // pointer to one histogram
		      int histo_size, // histo_size
		      int window_size, // size of the window_size in res value
		      int* start_interval_index,
		      int* end_interval_index,
		      int interval_lines);
void crosscorrelation_probability(int* clu1, // cells of interest
				  int* clu2,
				  int cell_pair_lines,
				  int* clu,  
				  int* res, 
				  int res_lines,
				  double* histo, // pointer to one histogram
				  int histo_size, // histo_size
				  int window_size, // size of the window_size in res value
				  int* start_interval_index,
				  int* end_interval_index,
				  int interval_lines);
SEXP crosscorrelation_cwrap(SEXP clu1_r, 
			    SEXP clu2_r,
			    SEXP cell_pair_lines_r,
			    SEXP clu_r,
			    SEXP res_r,
			    SEXP res_lines_r,
			    SEXP histo_size_r, // histo_size in bins
			    SEXP window_size_r, // size of the window_size in res value
			    SEXP start_interval_index_r,
			    SEXP end_interval_index_r,
			    SEXP interval_lines_r,
			    SEXP probability_r); // flag: 0 = count, 1 = probability
void cross_correlation_events(int* cells,
			      int cell_lines,
			      int* clu,  
			      int* res, 
			      int res_lines,
			      int* events,
			      int events_lines,
			      double* histo, // pointer to all histograms
			      int histo_size, // one histo size
			      int window_size, // size of the window_size in res value
			      int* start_interval_index, 
			      int* end_interval_index, 
			      int interval_lines);
SEXP crosscorrelationEvents_cwrap(SEXP cell_list_r, 
				  SEXP cell_list_lines_r,
				  SEXP clu_r,
				  SEXP res_r,
				  SEXP res_lines_r,
				  SEXP histo_size_r, // histo_size in bins
				  SEXP window_size_r, // size of the window_size in res value
				  SEXP start_interval_index_r,
				  SEXP end_interval_index_r,
				  SEXP interval_lines_r,
				  SEXP probability_r, // flag: 0 = count, 1 = probability
				  SEXP events_r,
				  SEXP events_lines_r) ;
void autocorrelation_one_cell(int clu_no, // cell of interest
			      int* clu,  
			      int* res, 
			      int res_lines,
			      int* histo, // pointer to one histogram
			      int histo_size, // histo_size
			      int window_size, // size of the window_size in res value
			      int* start_interval_index,
			      int* end_interval_index,
			      int interval_lines);
void autocorrelation_one_cell_probability(int clu_no, // cell of interest
					  int* clu,  
					  int* res, 
					  int res_lines,
					  double* histo, // pointer to one histogram
					  int histo_size, // histo_size
					  int window_size, // size of the window_size in res value
					  int* start_interval_index,
					  int* end_interval_index,
					  int interval_lines);
void autocorrelation(int* cell_list, // cells of interest
		     int cell_list_lines,
		     int* clu,  
		     int* res, 
		     int res_lines,
		     int* histo, // pointer to one histogram
		     int histo_size, // histo_size
		     int window_size, // size of the window_size in res value
		     int* start_interval_index,
		     int* end_interval_index,
		     int interval_lines);
void autocorrelation_probability(int* cell_list, // cells of interest
				 int cell_list_lines,
				 int* clu,  
				 int* res, 
				 int res_lines,
				 double* histo, // pointer to one histogram
				 int histo_size, // histo_size
				 int window_size, // size of the window_size in res value
				 int* start_interval_index,
				 int* end_interval_index,
				 int interval_lines);
SEXP autocorrelation_cwrap(SEXP cell_list_r, 
			   SEXP cell_list_lines_r,
			   SEXP clu_r,
			   SEXP res_r,
			   SEXP res_lines_r,
			   SEXP histo_size_r, // histo_size in bins
			   SEXP window_size_r, // size of the window_size in res value
			   SEXP start_interval_index_r,
			   SEXP end_interval_index_r,
			   SEXP interval_lines_r,
			   SEXP probability_r); // flag: 0 = count, 1 = probability


// intervals.c
void res_index_for_intervals(int* interval_lines, // value passed by reference in c
			     int* start, 
			     int* end, 
			     int res_lines, 
			     int* res, 
			     int* start_interval_index,
			     int* end_interval_index);
SEXP resIndexForIntervals_cwrap(SEXP interval_lines_r,
				SEXP start_r, 
				SEXP end_r, 
				SEXP res_lines_r, 
				SEXP res_r);

SEXP resWithinIntervals(SEXP interval_lines_r,
			SEXP start_r, 
			SEXP end_r, 
			SEXP res_lines_r, 
			SEXP res_r);



// firing_rate.c
void meanFiringRate(int* cells, // cells of interest
		    int cell_lines,
		    int* clu,  
		    int* res, 
		    int res_lines,
		    int* start_interval,
		    int* end_interval,
		    int* start_interval_index,
		    int* end_interval_index,
		    int interval_lines,
		    int sampling_rate,
		    double* rate);
SEXP meanFiringRate_cwrap(SEXP cell_list_r, 
			  SEXP cell_list_lines_r,
			  SEXP clu_r,
			  SEXP res_r,
			  SEXP res_lines_r,
			  SEXP start_interval_r,
			  SEXP end_interval_r,
			  SEXP start_interval_index_r,
			  SEXP end_interval_index_r,
			  SEXP interval_lines_r,
			  SEXP sampling_rate_r);

