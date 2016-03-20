#include <R.h>
#include <Rdefines.h>



/****************************************************************
Structures to read one data file                                     
****************************************************************/
struct data_file_si // data_file_short_integer
{
  char *file_name;
  int file_descriptor;
  int num_channels; 
  off_t file_size; 	                // length of the file in bytes
  long num_samples_in_file; 	// file_length/byte_per_sample
  short int* data_block;      // pointer to store the data from file
  int num_samples_in_complete_block; // number of samples in the complete blocks
  int block_size; // in bytes
};
/*************************************************************
Structure to read a group of dat file
**************************************************************/
struct group_data_file_si
{
  int num_files;
  int num_channels;
  struct data_file_si *file_group; // to hold the group of file
  long* resofs; // start sample of each file
  long num_samples_all_files;
  int last_file_index;
};

// dat_files.c
int init_data_file_si(struct data_file_si* df,const char *file_name,int num_channels);
int clean_data_file_si(struct data_file_si* df);
int data_file_si_load_block(struct data_file_si* df, long int start_index, long int size);
int data_file_si_get_data_one_channel(struct data_file_si* df, int channel_no, short int* one_channel, long int start_index, long int end_index);
int data_file_si_get_data_all_channels(struct data_file_si* df, short int* data, long int start_index, long int end_index);

int init_group_data_file_si(struct group_data_file_si* gdf, char** file_names,int num_files,int num_channels);
int group_data_file_si_get_data_one_channel(struct group_data_file_si* gf,int channel_no, short int* one_channel, long int start_index, long int end_index);
int clean_group_data_file_si(struct group_data_file_si* gdf);
SEXP group_data_file_si_get_one_channel_cwrap(SEXP file_names_r, SEXP num_channels_r, SEXP channel_no_r, SEXP start_index_r, SEXP end_index_r);





// math.c
double sum_double(int num_data, double* data, double invalid);
int find_max(int num_data, int* data);
double find_max_double_index(int num_data,double* data, int* index);
double find_max_double(int num_data,double* data);
double radian_to_degree(double radian);
double degree_to_radian(double degree);
void gaussian_kernel(double* kernel,int size, double standard_deviation);
void gaussian_kernel_2d(double* kernel,int x_size,int y_size,double standard_deviation);
void smooth_double_gaussian(double* array, int array_size, double smooth, double invalid);
SEXP smooth_double_gaussian_cwrap(SEXP array_r, SEXP array_size_r, SEXP sd_r, SEXP invalid_r);
void smooth_double_gaussian_degrees(double* array, int array_size,double smooth, double invalid);
SEXP smooth_double_gaussian_degrees_cwrap(SEXP array_r, SEXP array_size_r, SEXP sd_r, SEXP invalid_r);
void smooth_double_gaussian_circular(double* array, int array_size, double smooth, double invalid);
SEXP smooth_double_gaussian_circular_cwrap(SEXP array_r, SEXP array_size_r, SEXP sd_r, SEXP invalid_r);
void set_array_to_value_int (int* array, int array_size, int value);
void set_array_to_value_double (double* array, int array_size, double value);
void smooth_double_gaussian_2d(double* array, int x_size,int y_size, double smooth, double invalid);
SEXP smooth_double_gaussian_degrees_cwrap(SEXP array_r, SEXP array_size_r, SEXP sd_r, SEXP invalid_r);
double correlation (double* x, double* y, int size, double invalid);
void get_x_and_y_bin_from_index(int x_bins,int y_bins,int index,int* x_bin,int* y_bin);
int get_index_from_x_and_y_bin(int x_bins,int y_bins,int x_bin,int y_bin);
double distance(double x1, double y1, double x2, double y2);
SEXP detect_ttl_ups_cwrap(SEXP data_r,SEXP n_r,SEXP threshold_r);
SEXP detect_ttl_downs_cwrap(SEXP data_r,SEXP n_r,SEXP threshold_r);
void circular_stats_rate_histogram(double* histo, int num_bins, double* mean_direction, double* mean_vector_length);
SEXP circular_stats_rate_histogram_cwrap(SEXP cells_r, SEXP cell_lines_r,SEXP all_histos_r, SEXP histo_size_r);


// ifr.c
SEXP ifr_from_spike_density(SEXP res_r,SEXP clu_r, SEXP res_lines_r, SEXP window_size_ms_r, SEXP kernel_sd_ms_r, SEXP spike_bin_ms_r, SEXP cell_list_r, SEXP number_cells_r, SEXP start_interval_r, SEXP end_interval_r, SEXP interval_lines_r, SEXP sampling_rate_r);
void  firing_rate_per_cells_time_windows(int target_cell, int* res, int* clu, int res_lines, int bin_size_res,  double kernel_sd_res, double* firing_rate_in_bins, double* time_of_bin,int firing_rate_in_bins_lines,int* num_valid_bins, int* start_interval,int* end_interval, int* start_interval_index, int* end_interval_index,int interval_lines, int res_sampling_rate,int spike_bin_ms);


// convolution.c
void convolution_fftw3(double* inA, int inA_size, double* inB, int inB_size, double* out);
SEXP convolution_fftw3_cwrap(SEXP inA_r,SEXP inA_size_r, SEXP inB_r, SEXP inB_size_r);



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
double information_score(double* fr_map, double* occ_map, int map_size);
SEXP information_score_cwrap(SEXP cells_r,SEXP cell_lines_r, SEXP all_rate_maps_r,SEXP occupancy_map_r, SEXP map_size_r);
double sparsity_score(double* fr_map, double* occ_map, int map_size);
SEXP sparsity_score_cwrap(SEXP cells_r,SEXP cell_lines_r,SEXP all_rate_maps_r, SEXP occupancy_map_r, SEXP map_size_r);
SEXP map_autocorrelation_cwrap(SEXP cells_r, SEXP cell_lines_r,SEXP maps_r, SEXP num_bins_x_r, SEXP num_bins_y_r, SEXP auto_num_bins_x_r, SEXP auto_num_bins_y_r,SEXP min_bins_for_autocorrelation_r);
void map_autocorrelation(double *one_place, double *one_auto, int x_bins_place_map, int y_bins_place_map, int x_bins_auto_map, int y_bins_auto_map, int min_for_correlation);
void map_rotate(double* map,int x_bins,int y_bins,double deg,double invalid);
void detect_one_field_with_field( double* map, int x_bins, int y_bins, int min_num_bins_fields,double threshold, double* mean_x_field, double* mean_y_field, double* max_radius_field, int* num_bins_field,  double invalid,  double* field);
void detect_one_field( double* map, int x_bins, int y_bins, int min_num_bins_fields, double threshold, double* mean_x_field, double* mean_y_field, double* max_radius_field, int* num_bins_field, double invalid);
double hux_heading(double delta_x, double delta_y);
double degree_to_radian(double degree);
SEXP gridness_score_cwrap(SEXP cells_r, SEXP cell_lines_r, SEXP auto_maps_r, SEXP auto_num_bins_x_r, SEXP auto_num_bins_y_r, SEXP pixels_per_bin_r, SEXP number_fields_to_detect_r, SEXP min_num_bins_per_field_r, SEXP field_threshold_r, SEXP invalid_r);
double gridness_score(double* one_auto_map, int auto_num_bins_x,int auto_num_bins_y, double pixels_per_bin, int number_fields_to_detect,int min_num_bins_per_field, double field_threshold, double invalid);
double grid_orientation(double *map, int x_bins, int y_bins, double pixels_per_bin, int num_fields_to_detect, int min_num_bins_fields, float threshold, double invalid);
SEXP grid_orientation_cwrap(SEXP cells_r, SEXP cell_lines_r,SEXP auto_maps_r, SEXP auto_num_bins_x_r, SEXP auto_num_bins_y_r, SEXP pixels_per_bin_r,SEXP number_fields_to_detect_r, SEXP min_num_bins_per_field_r, SEXP field_threshold_r, SEXP invalid_r);
SEXP grid_spacing_cwrap(SEXP cells_r, SEXP cell_lines_r, SEXP auto_maps_r, SEXP auto_num_bins_x_r, SEXP auto_num_bins_y_r, SEXP pixels_per_bin_r, SEXP number_fields_to_detect_r, SEXP min_num_bins_per_field_r, SEXP field_threshold_r, SEXP invalid_r);
double grid_spacing(double *map, int x_bins, int y_bins, double pixels_per_bin, int num_fields_to_detect, int min_num_bins_fields, float threshold, double invalid);
int identify_border_pixels_in_occupancy_map(double* occ_map, int num_bins_x, int num_bins_y,int* border_map, int* border_x, int* border_y, int* num_bins_border);
int find_border_starting_point(double* occ_map, int num_bins_x, int num_bins_y,int*border_map,int*border_x,int* border_y,int* num_bins_border);
int find_an_adjacent_border_pixel(double* occ_map, int num_bins_x, int num_bins_y,int*border_map,int*border_x,int* border_y,int* num_bins_border);
void border_score_rectangular_environment(int* cells, int cell_lines, int num_bins_x, int num_bins_y, double* occ_map, double* maps, double percent_threshold_field, int min_bins_in_field, double* border_score, double* cm, double* dm, int* num_fields_detected);
SEXP border_score_rectangular_environment_cwrap(SEXP cells_r, SEXP cell_lines_r, SEXP num_bins_x_r,SEXP num_bins_y_r,SEXP occ_map_r,SEXP maps_r,SEXP percent_threshold_field_r,	SEXP min_bins_in_field_r);
void spike_head_direction(double *hd_whl, int whl_lines, int *res,  int res_lines,  double *hd_spike,  int res_samples_per_whl_sample, int* start_interval,  int* end_interval,  int interval_lines);
SEXP spike_head_direction_cwrap(SEXP hd_whl_r, SEXP whl_lines_r, SEXP res_r, SEXP res_lines_r, SEXP res_samples_per_whl_sample_r, SEXP start_interval_r, SEXP end_interval_r, SEXP interval_lines_r);
void occupancy_histogram(int x_bins,double pixels_per_bin_x, double *x_whl, int whl_lines, double *map, double ms_per_sample,int *start_interval, int *end_interval, int interval_lines, int res_samples_per_whl_sample, int num_repetitions);
SEXP occupancy_histogram_cwrap(SEXP x_bins_r, SEXP pixels_per_bin_x_r, SEXP x_whl_r, SEXP whl_lines_r, SEXP ms_per_sample_r, SEXP start_interval_r, SEXP end_interval_r, SEXP interval_lines_r, SEXP res_samples_per_whl_sample_r,SEXP n_repetitions_r);
void create_place_field_linear( int x_bins,double pixels_per_bin, double *x_spike, int *clu,int res_lines,int target_cell,double *occupancy_map,double *place_field,int num_repetitions);
SEXP firing_rate_histo_cwrap(SEXP num_bins_x_r, SEXP pixels_per_bin_x_r, SEXP x_spike_r, SEXP clu_r, SEXP res_lines_r, SEXP cells_r, SEXP num_cells_r, SEXP occ_histo_r, SEXP smooth_histo_sd_r, SEXP repetitions_r,SEXP circular_smooth_r);
SEXP detect_and_remove_field_cwrap(SEXP cells_r, SEXP cell_lines_r, SEXP map_r, SEXP x_bins_r, SEXP y_bins_r,SEXP num_fields_to_detect_r, SEXP min_num_bins_fields_r, SEXP threshold_r, SEXP invalid_r);
double detect_and_remove_field(double *map, int x_bins, int y_bins, int num_fields_to_detect,  int min_num_bins_fields,  float threshold, double invalid);
SEXP autocorrelation_doughnut_cwrap(SEXP cells_r,SEXP cell_lines_r, SEXP maps_r,SEXP x_bins_r,SEXP y_bins_r,SEXP num_fields_to_detect_r,SEXP min_num_bins_fields_r,SEXP threshold_r,SEXP px_per_bin_r,SEXP invalid_r);
void autocorrelation_doughnut(double* one_auto_map, int auto_num_bins_x,int auto_num_bins_y,  double pixels_per_bin,int number_fields_to_detect, int min_num_bins_per_field,double field_threshold, double invalid);
SEXP autocorrelation_doughnut_rotate_cwrap(SEXP cells_r,SEXP cell_lines_r, SEXP maps_r, SEXP x_bins_r, SEXP y_bins_r, SEXP num_fields_to_detect_r, SEXP min_num_bins_fields_r,SEXP threshold_r,SEXP px_per_bin_r,SEXP rotations_r,SEXP degree_r,SEXP invalid_r);
void autocorrelation_doughnut_rotate(double* one_auto_map, int auto_num_bins_x, int auto_num_bins_y, double* rotated_maps,double pixels_per_bin,int number_fields_to_detect, int min_num_bins_per_field, double field_threshold, int rotations, double degree, double invalid);
void spike_position_1d(double *x_whl,int whl_lines,int *res,int res_lines, double *x_spike, int res_samples_per_whl_sample,int* start_interval,
int* end_interval, int interval_lines);
SEXP spike_position_1d_cwrap(SEXP x_whl_r, SEXP whl_lines_r, SEXP res_r, SEXP res_lines_r,SEXP res_samples_per_whl_sample_r,SEXP start_interval_r,
SEXP end_interval_r, SEXP interval_lines_r);


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

