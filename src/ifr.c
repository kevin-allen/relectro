#include <math.h>
#include "relectro.h"

SEXP ifr_from_spike_density(SEXP res_r,SEXP clu_r, SEXP res_lines_r, 
			    SEXP window_size_ms_r, SEXP kernel_sd_ms_r, SEXP spike_bin_ms_r,
			    SEXP cell_list_r, SEXP cell_lines_r,
			    SEXP start_interval_r, SEXP end_interval_r, SEXP interval_lines_r,
			    SEXP sampling_rate_r)
{
  int window_size_res= (int)(REAL(window_size_ms_r)[0]*REAL(sampling_rate_r)[0]/1000.0);
  double kernel_sd_res = REAL(kernel_sd_ms_r)[0]*REAL(sampling_rate_r)[0]/1000.0;
  int interval_lines=INTEGER_VALUE(interval_lines_r);
  int res_lines=INTEGER_VALUE(res_lines_r);
  int cell_lines=INTEGER_VALUE(cell_lines_r);
  int* res = INTEGER_POINTER(res_r);
  int* cell_list = INTEGER_POINTER(cell_list_r);
  int number_windows=0;
  int num_valid_bins=0;
  
  // set maximum to 
  int max_tp, max_res, max_int;
  max_res=find_max(res_lines,res);
  max_int=find_max(interval_lines,INTEGER_POINTER(end_interval_r));
  if(max_res>max_int)
    {
     max_tp=max_res;
   }
  else{
    max_tp=max_int;
  }
  
  
  number_windows=(max_tp/window_size_res);
  if(max_tp%window_size_res!=0)
    {number_windows++;}
  if (number_windows <1)
    {
      return 0;
    }

  int* start_interval_index = (int*) malloc (interval_lines*sizeof(int));
  int* end_interval_index =  (int*) malloc (interval_lines*sizeof(int));
  res_index_for_intervals(&interval_lines,
			       INTEGER_POINTER(start_interval_r),
			       INTEGER_POINTER(end_interval_r), 
			       res_lines, 
			       res,
			       start_interval_index,
			       end_interval_index,0);

 number_windows=(max_tp/window_size_res)+1;

 //Rprintf("number_windows: %d\n",number_windows);
 //Rprintf("max_tp: %d\n",max_tp);
 
 double* instantaneous_fr= (double*) malloc(number_windows*cell_lines*sizeof(double)); // cells x windows
 double* time_of_bin= (double*) malloc(number_windows*sizeof(double)); // windows
 double* ptr_ifr;
 int target_cell;
 for(int i = 0; i < cell_lines; i++)
 {
  target_cell=cell_list[i];
  ptr_ifr=instantaneous_fr+(number_windows*i);
  firing_rate_per_cells_time_windows(target_cell,
                                     res, 
                                     INTEGER_POINTER(clu_r), 
                                     res_lines,
                                     window_size_res, 
                                     kernel_sd_res, 
                                     ptr_ifr,  // allocated for the longest size possible
                                     time_of_bin, // allocated for the longest size possible
                                     number_windows, // longest size possible with all cells
                                     &num_valid_bins, // actual size, single integer that will be modified
                                     INTEGER_POINTER(start_interval_r),
                                     INTEGER_POINTER(end_interval_r),
                                     start_interval_index,
                                     end_interval_index,
                                     interval_lines,
                                     REAL(sampling_rate_r)[0],
                                     REAL(spike_bin_ms_r)[0]);
 }
 
 
 SEXP ans;
 double* rans;
 PROTECT(ans = allocMatrix(REALSXP, cell_lines+1, num_valid_bins)); // x y, matrix[x,y], x is a Row in R
 rans = REAL(ans);
 
 for(int j = 0; j < num_valid_bins; j++){
   rans[(cell_lines+1)*j]=time_of_bin[j];
 }
 
 for(int i = 1 ; i < cell_lines+1; i++){
   ptr_ifr=instantaneous_fr+(number_windows*(i-1)); // pointer to cell 0
 for(int j = 0; j < num_valid_bins; j++){
   rans[(cell_lines+1)*j+i]=ptr_ifr[j];
  }
 }
 

 free(instantaneous_fr);
 free(time_of_bin);
 free(start_interval_index);
 free(end_interval_index);
 UNPROTECT(1);
 return(ans);
}

void  firing_rate_per_cells_time_windows(int target_cell,
                                         int* res, 
					 int* clu, 
					 int res_lines, 
					 int window_size_res, 
					 double kernel_sd_res, 
					 double* firing_rate_in_bins,  // allocated for the longest size possible
					 double* time_of_bin, // allocated for the longest size possible
					 int firing_rate_in_bins_lines, // longest size possible
					 int* num_valid_bins, // actual size, single integer that will be modified
					 int* start_interval,
					 int* end_interval,
					 int* start_interval_index,
					 int* end_interval_index,
					 int interval_lines,
					 int res_sampling_rate,
					 int spike_bin_ms)
{
  /*
    To get an estimate of local spike density of
    each cell i as a function of time bin t,
    spike trains were convolved with Gaussian kernels 
    and binned at x ms
    
    Similar to what Durstewitz seems to do
    
    calculate local spike density at each ms and then integrate it for 
  */

  int res_per_data_point_to_fft=spike_bin_ms*res_sampling_rate/1000;
  int win=0;
  int start_bin, end_bin;
  int num_bins=firing_rate_in_bins_lines;
  double* pt;
  double* pt_time;
  double* local_spike_density; // each time point get local spike density from convol.
  double* kernel;
  double* spike_array;
  int lsp_start;
  int lsp_end;
  int num_standard_deviations_in_kernel=5;
  set_array_to_value_double (firing_rate_in_bins,firing_rate_in_bins_lines,0);
  set_array_to_value_double (time_of_bin,firing_rate_in_bins_lines,-1);
 
 
 
 int max_tp, max_res, max_int;
 max_res=find_max(res_lines,res);
 max_int=find_max(interval_lines,end_interval);
 if(max_res>max_int)
 {
   max_tp=max_res;
 }
 else{
   max_tp=max_int;
 }
 
  // create a gaussian kernel for convolution
  int kernel_size=((int)((kernel_sd_res/res_per_data_point_to_fft)*num_standard_deviations_in_kernel*2))+1;
  if (kernel_size%2!=1) // should be an odd number
    {
      kernel_size++;
    }
  kernel=(double*) malloc(kernel_size*sizeof(double));
  gaussian_kernel(kernel,
		  kernel_size,
		  kernel_sd_res/res_per_data_point_to_fft);
  
  
  // power of 2 are done in convolution function
  int size_for_fft=max_tp/res_per_data_point_to_fft+1;
  
  if(kernel_size>size_for_fft)
  {
    Rprintf("firing_rate_per_cells_time_windows: kernel_size>size_for_fft\n");
    return;
  }
    
  local_spike_density = (double*) malloc(size_for_fft*sizeof(double)); // array for all samples
  spike_array =  (double*) malloc(size_for_fft*sizeof(double));
  
  
  // get the spikes in local spike density array
  set_array_to_value_double(spike_array,size_for_fft,0);
  for (int j = 0; j < res_lines;j++)
  {
	  if(clu[j]==target_cell&&res[j]!=-1)
	    spike_array[res[j]/res_per_data_point_to_fft]=1;
	}
  
  /*
  Rprintf("target_cell:%d\n",target_cell);
  Rprintf("res_per_data_point_to_fft:%d\n",res_per_data_point_to_fft);
  Rprintf("firing_rate_in_bins_lines:%d\n",firing_rate_in_bins_lines);
  Rprintf("kernel_size:%d\n",kernel_size);
  Rprintf("size_for_fft:%d\n",size_for_fft);
  */
  
  
  
  // do the convolution of local_spike_density with padded_kernel    
  convolution_fftw3(spike_array,
	                  size_for_fft,
	                  kernel,
	                  kernel_size,
	                  local_spike_density); // one data point each 20 res value
      
  // get pointers for the data of this cell
  
  
  win=0;
  // loop for each interval
  for (int j = 0; j < interval_lines;j++)
 	{
	  start_bin=start_interval[j];
	  end_bin=start_bin+window_size_res;
	  while(end_bin<end_interval[j]) // all the bin needs to be in interval
	    {
	      // integrate for each bin the local spike density
	      lsp_start=start_bin/res_per_data_point_to_fft;
	      lsp_end=end_bin/res_per_data_point_to_fft;
	      for(int k = lsp_start;k<lsp_end;k++)
	  	  {
	        firing_rate_in_bins[win]+=local_spike_density[k];
		    } 
	      firing_rate_in_bins[win]=firing_rate_in_bins[win]*(res_sampling_rate/window_size_res); // Hz
	      time_of_bin[win]=(double)(start_bin+(end_bin-start_bin)/2)/res_sampling_rate; //Sec
	      win++;
	      start_bin=end_bin;
	      end_bin=start_bin+window_size_res;
	    }
	}
 // Rprintf("target_cell: %d\n",target_cell);
//  Rprintf("win: %d\n",win);
//  Rprintf("firing_rate_in_bins_lines: %d\n",firing_rate_in_bins_lines);
//  Rprintf("max fr: %lf\n",find_max_double(win,firing_rate_in_bins));
  
  *num_valid_bins=win;
  free(local_spike_density);
  free(kernel);
  free(spike_array);
  return;
}
