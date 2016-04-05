#include "relectro.h"


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
		    double* rate)
{
  set_array_to_value_double(rate,cell_lines,0);
  
  // count the spikes in intervals
  for (int j = 0; j < interval_lines ; j++)
    for (int i = 0; i < cell_lines; i++)
    {
      for (int x = start_interval_index[j]; x <= end_interval_index[j]; x++)
        if (clu[x] == cells[i])
          rate[i]++;
    }
  int time_res=0;
  double time_sec=0;
  for (int i = 0; i < interval_lines ; i++)
    time_res = time_res + (end_interval[i]-start_interval[i]);
  time_sec = (double)time_res/(double)sampling_rate; // in hz
  for (int i = 0; i < cell_lines; i++){
    rate[i]=rate[i]/time_sec;
  }
  return;
}


// wrapper to call firing rate function from .Call()
// same name as c function with _cwrap suffix
// all object types are SEXP
// argument names have the suffix _r
// should deal with count and probability
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
			  SEXP sampling_rate_r)
{ 
  // transform the SEXP object in their correct c types
  // create the list of variable of correct c types
  int* cell_list;
  int cell_list_lines;
  int* clu;
  int* res;
  int res_lines;
  int* start_interval;
  int* end_interval;
  int* start_interval_index;
  int* end_interval_index;
  int interval_lines;
  int sampling_rate;
  
  /* // coersion */
  cell_list=INTEGER_POINTER(cell_list_r);
  cell_list_lines=INTEGER_VALUE(cell_list_lines_r); //get an integer
  clu=INTEGER_POINTER(clu_r);
  res=INTEGER_POINTER(res_r);
  res_lines=INTEGER_VALUE(res_lines_r);
  start_interval=INTEGER_POINTER(start_interval_r);
  end_interval=INTEGER_POINTER(end_interval_r);
  start_interval_index=INTEGER_POINTER(start_interval_index_r);
  end_interval_index=INTEGER_POINTER(end_interval_index_r);
  interval_lines=INTEGER_VALUE(interval_lines_r);
  sampling_rate=INTEGER_VALUE(sampling_rate_r);
  
  // prepare a SEXP object with the data from the histogram
  SEXP out = PROTECT(allocVector(REALSXP, cell_list_lines));
  // allocate memory for histo
  double* rate = REAL(out);
  // call the c function from electrophys library
  meanFiringRate(cell_list, // cells of interest
		 cell_list_lines,
		 clu,
		 res,
		 res_lines,
		 start_interval,
		 end_interval,
		 start_interval_index,
		 end_interval_index,
		 interval_lines,
		 sampling_rate,
		 rate);
  //free memory for rate
  UNPROTECT(1);
  return(out);
}

