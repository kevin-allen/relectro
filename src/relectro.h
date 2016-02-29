#include <R.h>
#include <Rdefines.h>

int find_max(int num_data, int* data);
int check_interval_chronology_between(int num_lines, 
				      int* start, 
				      int* end);

void set_res_outside_interval_to_minus_one(int interval_lines,
					   int* start_interval_index,
					   int* end_interval_index,
					   int res_lines,
					   int* res);
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

