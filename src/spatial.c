#include <math.h>
#include "relectro.h"


SEXP angular_speed_from_hd_cwrap(SEXP hd_r,
			  SEXP whl_lines_r,
			  SEXP look_back_max_r,
			  SEXP look_ahead_max_r,
			  SEXP res_sampling_rate_r, 
			  SEXP res_samples_per_whl_sample_r)

{
  int whl_lines=INTEGER_VALUE(whl_lines_r);
  SEXP out = PROTECT(allocVector(REALSXP, whl_lines));
  double* speed = (double*)malloc(whl_lines*sizeof(double));

  // call c function
  angular_speed_from_hd(REAL(hd_r),
  		 speed,
  		 whl_lines,
  		 INTEGER_VALUE(look_back_max_r),
  		 INTEGER_VALUE(look_ahead_max_r),
  		 INTEGER_VALUE(res_sampling_rate_r),
  		 INTEGER_VALUE(res_samples_per_whl_sample_r));

  for(int i=0;i< whl_lines;i++)
    REAL(out)[i]=speed[i];
  free(speed);
  UNPROTECT(1);
  return(out);
}

void angular_speed_from_hd(double* hd, double* angular_speed, int whl_lines, int look_back_max, int look_ahead_max, int res_sampling_rate, int res_samples_per_whl_sample)
{
  double samples_per_sec=(double)res_sampling_rate/(double)res_samples_per_whl_sample;
  int after_distance = 0;
  int before_distance = 0;
  int before_index=0;
  int after_index=0;
  double diff_deg;
  for (int i = 0; i < whl_lines; i++)
   {
     if ( hd[i] == -1.0)
       {
	 angular_speed[i]=-1.0;
       }
     else
       {
	 after_index=i;
	 after_distance=0;
	 before_index=i;
	 before_distance=0;
	 // look after except if the last value in whl
	 if (i < whl_lines-1)
	   {
	     // look if valid ok after
	     do
	       {
		 after_index++;
		 after_distance++;
	       } while ((hd[after_index] == -1.0) && 
			(after_distance < look_ahead_max) &&
			(after_index < whl_lines));
	   }
	 // look before except if the first value in whl
	 if ( i > 0)
	   {
	     // look if valid ok before
	     do
	       {
		 before_index--;
		 before_distance++;
	       } while ((hd[before_index] == -1.0) && 
			(before_distance < look_back_max) &&
			(before_index >= 0));
	   }

	  if(hd[before_index] != -1.0 && hd[after_index] != -1.0)
	   {
	     //get the difference in degrees
	     if(hd[after_index]>hd[before_index])
	       {
		 diff_deg=hd[after_index]-hd[before_index];
		 if(diff_deg>180)
		   diff_deg=360-diff_deg;
	       }
	     else
	       {
		 diff_deg=hd[before_index]-hd[after_index];
		 if(diff_deg>180)
		   diff_deg=360-diff_deg;
	       }
	     diff_deg=diff_deg/(after_distance+before_distance); // take into account time before and after, so we get deg per sample
	     angular_speed[i] = diff_deg*samples_per_sec;
	   }
	  else
	    {
	      angular_speed[i] = -1.0;
	    }
       }
   }
}

void speed_from_whl(double* x_whl,
		    double* y_whl,
		    double* speed,
		    int whl_lines,
		    int look_back_max,
		    int look_ahead_max,
		    double px_per_cm,
		    int res_sampling_rate, // ex 20 000
		    int res_samples_per_whl_sample) // ex 512
{
  double samples_per_sec=(double)res_sampling_rate/(double)res_samples_per_whl_sample;
  int after_distance = 0;
  int before_distance = 0;
  int before_index=0;
  int after_index=0;
  double diff_x;
  double diff_y;
  double distance_px;
  double distance_cm;
  for (int i = 0; i < whl_lines; i++)
   {
     if ( x_whl[i] == -1.0 || y_whl[i] == -1.0)
       {
	 speed[i]=-1.0;
       }
     else
       {
	 after_index=i;
	 after_distance=0;
	 before_index=i;
	 before_distance=0;

	 // look after except if the last value in whl
	 if (i < whl_lines-1)
	   {
	     // look if valid ok after
	     do
	       {
		 after_index++;
		 after_distance++;
	       } while ((x_whl[after_index] == -1.0) && 
			(after_distance < look_ahead_max) &&
			(after_index < whl_lines));
	   }
	 // look before except if the first value in whl
	 if ( i > 0)
	   {
	     // look if valid ok before
	     do
	       {
		 before_index--;
		 before_distance++;
	       } while ((x_whl[before_index] == -1.0) && 
			(before_distance < look_back_max) &&
			(before_index >= 0));
	   }
	  if(x_whl[before_index] != -1.0 && x_whl[after_index] != -1.0)
	   {
	     // get the distance between current and after
	     diff_x = sqrt(pow((x_whl[before_index] - x_whl[after_index]),2)); // distance px x
	     diff_y = sqrt(pow((y_whl[before_index] - y_whl[after_index]),2)); // distance px y
	     distance_px = diff_x*diff_x + diff_y*diff_y; // pythagoras
	     distance_px = sqrt(distance_px)/ (after_distance+before_distance); //  pythagoras + ahead or back
	     distance_cm = distance_px/px_per_cm; // in cm
	     speed[i] = distance_cm*samples_per_sec;
	   }
	  else
	    {
	      speed[i] = -1.0;
	    }
       }
   }
}

SEXP speed_from_whl_cwrap(SEXP x_whl_r,
			  SEXP y_whl_r,
			  SEXP whl_lines_r,
			  SEXP look_back_max_r,
			  SEXP look_ahead_max_r,
			  SEXP px_per_cm_r,
			  SEXP res_sampling_rate_r, 
			  SEXP res_samples_per_whl_sample_r)

{
  int whl_lines=INTEGER_VALUE(whl_lines_r);
  SEXP out = PROTECT(allocVector(REALSXP, whl_lines));
  double* speed = (double*)malloc(whl_lines*sizeof(double));

  // call c function
  speed_from_whl(REAL(x_whl_r),
  		 REAL(y_whl_r),
  		 speed,
  		 whl_lines,
  		 INTEGER_VALUE(look_back_max_r),
  		 INTEGER_VALUE(look_ahead_max_r),
  		 REAL(px_per_cm_r)[0],
  		 INTEGER_VALUE(res_sampling_rate_r),
  		 INTEGER_VALUE(res_samples_per_whl_sample_r));

  for(int i=0;i< whl_lines;i++)
    REAL(out)[i]=speed[i];
  
  free(speed);
  UNPROTECT(1);
  return(out);
}

SEXP speed_at_res_values_cwrap(SEXP speed_r,
			       SEXP whl_lines_r,
			       SEXP res_r,
			       SEXP res_lines_r,
			       SEXP res_samples_per_whl_sample_r)
{
  int whl_lines=INTEGER_VALUE(whl_lines_r);
  int res_lines=INTEGER_VALUE(res_lines_r);
  SEXP out = PROTECT(allocVector(REALSXP, res_lines));
  double* speed_at_res = (double*)malloc(res_lines*sizeof(double));

  spike_position_no_interval(REAL(speed_r),
			     whl_lines,
			     INTEGER_POINTER(res_r),
			     res_lines,
			     speed_at_res,
			     INTEGER_VALUE(res_samples_per_whl_sample_r));

  for(int i=0;i< res_lines;i++)
    REAL(out)[i]=speed_at_res[i];
  
  free(speed_at_res);
  UNPROTECT(1);
  return(out);



}

void spike_position_no_interval(double *x_whl,
				int whl_lines,
				int *res,
				int res_lines,
				double *x_spike,
				int res_samples_per_whl_sample)
{ //Function to get the linear position of each spike
  //Can be used to get the speed at given res values
  int whl_index;
  int residual;
  for (int i = 0; i < res_lines; i++)
    {
      if (res[i] != -1.0)
	{
	  whl_index = (res[i] / res_samples_per_whl_sample) - 1;
	  if (whl_index >=  whl_lines || (whl_index < 0 && whl_index != -1))
	    {
	      printf("whl file is too small for some res values:%d\n",res[i]);
	    }
	  residual = res[i] % res_samples_per_whl_sample;
	  // spikes before the first frame capture
	  if (whl_index == -1)
	    {
	      x_spike[i] = x_whl[0];
	    }
	  
	  // spikes between the first and last frame capture
	  if ((whl_index < whl_lines -1) && (whl_index != -1))
	    {
	      x_spike[i] =
		x_whl[whl_index] +
		((x_whl[whl_index+1]- x_whl[whl_index])
		 /res_samples_per_whl_sample * residual);
	    }
	  
	  // if spikes after the last frame capture
	  if (whl_index == whl_lines -1)
	    {
	      x_spike[i] = x_whl[whl_index];
	    }
	  
	  
	  // cope with the -1 in whl
	  if (x_whl[whl_index] == -1.0 || x_whl[whl_index+1]==-1.0)
	    {
	      x_spike[i] = -1.0;
	    }
	}
      else // the res value was set at -1
	{
	  x_spike[i] = -1.0;
	}
    }
  return;
}



SEXP speed_intervals_cwrap(SEXP speed_r, SEXP whl_lines_r, SEXP res_samples_per_whl_sample_r, SEXP min_speed_r, SEXP max_speed_r)
{
  int whl_lines=INTEGER_VALUE(whl_lines_r);
  int interval_lines;
  interval_lines=speed_intervals_count(REAL(speed_r), 
				       whl_lines, 
				       INTEGER_VALUE(res_samples_per_whl_sample_r),
				       REAL(min_speed_r)[0],
				       REAL(max_speed_r)[0]);
  if(interval_lines<1)
    {
      printf("speed_intervals_cwrap: number of intervals < 1");
      return(R_NilValue);
    }
  int* start = (int*)malloc(interval_lines*sizeof(int));
  int* end = (int*)malloc(interval_lines*sizeof(int));
  SEXP out = PROTECT(allocMatrix(REALSXP,2,interval_lines));
  speed_intervals(REAL(speed_r), 
		  whl_lines, 
		  INTEGER_VALUE(res_samples_per_whl_sample_r),
		  REAL(min_speed_r)[0],
		  REAL(max_speed_r)[0],
		  start,
		  end);
  double* rans;
  rans = REAL(out);
  // copy the results in a SEXP
  for(int i = 0; i < interval_lines; i++) {
    rans[i*2] = start[i];
    rans[i*2+1]=end[i];
  }
  free(start);
  free(end);
  UNPROTECT(1);
  return(out);
}



int speed_intervals_count(double* speed, int whl_lines, int res_samples_per_whl_sample,double min_speed,double max_speed)
{ // returns how many intervals will be found
  int minus_one=0;
  int was_in=0;
  int max_minus_one_before_out=10;
  int end_int=-1;
  int start_int=-1;
  int jj=0;
  int num_intervals=0;
  for (int i = 0; i < whl_lines; i++)
    {
      if (speed[i]==-1)
	{
	  minus_one++;
	  if (was_in==1 && minus_one > max_minus_one_before_out)
	    { // to exit of interval when we lose track of animal
	      end_int=(i*res_samples_per_whl_sample)+res_samples_per_whl_sample;
	      was_in=0;
	      jj=0;
	      num_intervals++;
	      //	      cout << start_int << " " << end_int << '\n';
	    }
	}
      if (speed[i]!=-1)
	{
	  minus_one=0;
	  if (was_in==1) // look for exit time and print
	    {
	      if (speed[i]<min_speed||speed[i]>max_speed)
		{
		  end_int=(i*res_samples_per_whl_sample)+res_samples_per_whl_sample;
		  jj=0;
		  //  cout << start_int << " " << end_int << '\n';
		  num_intervals++;
		}
	    }
	  if (was_in==0) // look for entry time in the speed interval
	    {
	      if(speed[i]>=min_speed&&speed[i]<=max_speed)
		{
		  start_int=(i*res_samples_per_whl_sample)+res_samples_per_whl_sample;
		  jj=1;
		}
	    }
	  was_in=jj;
	}
      
      // to print exit time if the animal is within speed interval and reach the end of file
      if (i==whl_lines-1 && was_in==1)
	{
	  end_int=(i*res_samples_per_whl_sample)+res_samples_per_whl_sample;
	  if (start_int != end_int)
	    {
	      jj=0;
	      num_intervals++;
	      //	      cout << start_int << " " << end_int << '\n';
	    }
	}

    }
  return(num_intervals);
}
void speed_intervals(double* speed, int whl_lines, int res_samples_per_whl_sample,double min_speed,double max_speed,int* start,int* end)
{ // returns the intervals
  int minus_one=0;
  int was_in=0;
  int max_minus_one_before_out=10;
  int end_int=-1;
  int start_int=-1;
  int jj=0;
  int num_intervals=0;
  for (int i = 0; i < whl_lines; i++)
    {
      if (speed[i]==-1)
	{
	  minus_one++;
	  if (was_in==1 && minus_one > max_minus_one_before_out)
	    { // to exit of interval when we lose track of animal
	      end_int=(i*res_samples_per_whl_sample)+res_samples_per_whl_sample;
	      was_in=0;
	      jj=0;
	      start[num_intervals]=start_int;
	      end[num_intervals]=end_int;
	      num_intervals++;
	      //	      cout << start_int << " " << end_int << '\n';
	    }
	}
      if (speed[i]!=-1)
	{
	  minus_one=0;
	  if (was_in==1) // look for exit time and print
	    {
	      if (speed[i]<min_speed||speed[i]>max_speed)
		{
		  end_int=(i*res_samples_per_whl_sample)+res_samples_per_whl_sample;
		  jj=0;
		  start[num_intervals]=start_int;
		  end[num_intervals]=end_int;
		  //  cout << start_int << " " << end_int << '\n';
		  num_intervals++;
		}
	    }
	  if (was_in==0) // look for entry time in the speed interval
	    {
	      if(speed[i]>=min_speed&&speed[i]<=max_speed)
		{
		  start_int=(i*res_samples_per_whl_sample)+res_samples_per_whl_sample;
		  jj=1;
		}
	    }
	  was_in=jj;
	}
      
      // to print exit time if the animal is within speed interval and reach the end of file
      if (i==whl_lines-1 && was_in==1)
	{
	  end_int=(i*res_samples_per_whl_sample)+res_samples_per_whl_sample;
	  if (start_int != end_int)
	    {
	      jj=0;
	      start[num_intervals]=start_int;
	      end[num_intervals]=end_int;
	      num_intervals++;
	      //	      cout << start_int << " " << end_int << '\n';
	    }
	}

    }
  return;
}









SEXP spike_position_cwrap(SEXP x_whl_r,
			  SEXP y_whl_r,
			  SEXP whl_lines_r,
			  SEXP res_r,
			  SEXP res_lines_r,
			  SEXP res_samples_per_whl_sample_r,
			  SEXP start_interval_r,
			  SEXP end_interval_r,
			  SEXP interval_lines_r)
{
  int whl_lines=INTEGER_VALUE(whl_lines_r);
  int res_lines=INTEGER_VALUE(res_lines_r);
  SEXP out = PROTECT(allocMatrix(REALSXP,2,res_lines));
  double* x_spike = (double*)malloc(res_lines*sizeof(double));
  double* y_spike = (double*)malloc(res_lines*sizeof(double));
  spike_position(REAL(x_whl_r),
		 REAL(y_whl_r),
		 whl_lines,
		 INTEGER_POINTER(res_r),
		 res_lines,
		 x_spike,
		 y_spike,
		 INTEGER_VALUE(res_samples_per_whl_sample_r),
		 INTEGER_POINTER(start_interval_r),
		 INTEGER_POINTER(end_interval_r),
		 INTEGER_VALUE(interval_lines_r));
  double* rans;
  rans = REAL(out);
  // copy the results in a SEXP
  for(int i = 0; i < res_lines; i++) {
    rans[i*2] = x_spike[i];
    rans[i*2+1]= y_spike[i];
  }
  free(x_spike);
  free(y_spike);
  UNPROTECT(1);
  return(out);

}



void spike_position(double *x_whl,
		    double *y_whl,
		    int whl_lines,
		    int *res,
		    int res_lines,
		    double *x_spike,
		    double *y_spike,
		    int res_samples_per_whl_sample,
		    int* start_interval,
		    int* end_interval,
		    int interval_lines)
{
  ///////////////////////////////////////////////////////////////////
  //	Function to get the position of each spike within intervals//
  ///////////////////////////////////////////////////////////////////

  // get the res index for the intervals
  int* start_interval_index;
  int* end_interval_index;
  start_interval_index = (int*)malloc(interval_lines*sizeof(int));
  end_interval_index = (int*)malloc(interval_lines*sizeof(int));

  res_index_for_intervals(&interval_lines,
			  start_interval,
			  end_interval, 
			  res_lines, 
			  res,
			  start_interval_index,
			  end_interval_index);
  // set all the spike position to -1
  set_array_to_value_double(x_spike,res_lines,-1.0);
  set_array_to_value_double(y_spike,res_lines,-1.0);

  //loop for each interval and set the position when valid//
  int whl_index;
  int residual;
    for (int k = 0; k < interval_lines; k++)
    {
      for (int i = start_interval_index[k]; i < end_interval_index[k]; i++)
	{
	  if (res[i] != -1.0)
	    {
	      whl_index = (res[i] / res_samples_per_whl_sample) - 1;
	      if (whl_index >=  whl_lines || (whl_index < 0 && whl_index != -1))
		{
		  printf("whl file is too small for some res value\n");
		}
	      residual = res[i] % res_samples_per_whl_sample;
	      // spikes before the first frame capture
	      if (whl_index == -1)
		{
		  x_spike[i] = x_whl[0];
		  y_spike[i] = y_whl[0];
		}
	      // spikes between the first and last frame capture
	      if ((whl_index < whl_lines -1) && (whl_index != -1))
		{
		  x_spike[i] =
		    x_whl[whl_index] +
		    ((x_whl[whl_index+1]- x_whl[whl_index])
		     /res_samples_per_whl_sample * residual);
		  
		  y_spike[i] =
		    y_whl[whl_index] +
		    ((y_whl[whl_index+1]- y_whl[whl_index])
		     /res_samples_per_whl_sample * residual);
		}
	      
	      // if spikes after the last frame capture
	      if (whl_index == whl_lines -1)
		{
		  x_spike[i] = x_whl[whl_index];
		  y_spike[i] = y_whl[whl_index];
		}
	      // cope with the -1 values from the whl file
	      if (x_whl[whl_index] == -1.0 || x_whl[whl_index+1]==-1.0)
		{
		  x_spike[i] = -1.0;
		  y_spike[i] = -1.0;
		}
	    }
	  else // the res value was set at -1
	    {
	      x_spike[i] = -1.0;
	      y_spike[i] = -1.0;
	    }
	}
    }
    free(start_interval_index);
    free(end_interval_index);
    return;
}

SEXP occupancy_map_cwrap(SEXP x_bins_r, // num of bins x
			 SEXP y_bins_r, // num of bins y
			 SEXP pixels_per_bin_x_r,
			 SEXP pixels_per_bin_y_r,
			 SEXP x_whl_r, 
			 SEXP y_whl_r,
			 SEXP whl_lines_r, // num lines in whl data
			 SEXP ms_per_sample_r, // ms per sample for whl data
			 SEXP start_interval_r, // starts of the intervals in res value
			 SEXP end_interval_r, // ends of intervals in res value
			 SEXP interval_lines_r, // number of intervals
			 SEXP res_samples_per_whl_sample_r)
{

  int x_bins=INTEGER_VALUE(x_bins_r);
  int y_bins=INTEGER_VALUE(y_bins_r);
  SEXP out = PROTECT(allocMatrix(REALSXP,y_bins,x_bins)); // nrow ncol, to return
  double* map = (double*)malloc(y_bins*x_bins*sizeof(double)); // to pass to c function
 
  occupancy_map(x_bins, // num of bins x
		y_bins, // num of bins y
		REAL(pixels_per_bin_x_r)[0],
		REAL(pixels_per_bin_y_r)[0],
		REAL(x_whl_r),
		REAL(y_whl_r),
		INTEGER_VALUE(whl_lines_r), // num lines in whl data
		map, // occupancy map
		REAL(ms_per_sample_r)[0], // ms per sample for whl data
		INTEGER_POINTER(start_interval_r), // starts of the intervals in res value
		INTEGER_POINTER(end_interval_r), // ends of intervals in res value
		INTEGER_VALUE(interval_lines_r), // number of intervals
		INTEGER_VALUE(res_samples_per_whl_sample_r));

  double* rans;
  rans = REAL(out);
  // copy the results in a SEXP
  for(int x = 0; x < x_bins; x++)
    for(int y = 0; y < y_bins; y++){
      rans[x*y_bins+y] = map[x*y_bins+y];
    }
  free(map);
  UNPROTECT(1);
  return(out);
}

void occupancy_map(int x_bins, // num of bins x
		   int y_bins, // num of bins y
		   double pixels_per_bin_x,
		   double pixels_per_bin_y,
		   double *x_whl, 
		   double *y_whl,
		   int whl_lines, // num lines in whl data
		   double *map, // occupancy map
		   double ms_per_sample, // ms per sample for whl data
		   int *start_interval, // starts of the intervals in res value
		   int *end_interval, // ends of intervals in res value
		   int interval_lines, // number of intervals
		   int res_samples_per_whl_sample)
{
  /**********************************************************************
     Create the occupancy map for 2d firing map
     This add  ms_per_sample for every valid whl value that is in 
     the intervals given valid whl values outside the intervals 
     will not be added to occupancy map.
     If a valid whl value is for a period outside intervals,
     the right proportion of time will be removed
    
     There is a minimum time that the animal has to 
     spend in a bin for it to be valid; see variable min_occ

     Bins that are on there own in space are removed.
  ************************************************************************/
  int min_occ = 20; // minimum time spent in a bin
    // variable to deal with the intervals
  int res_value_for_whl_sample; // time in res associated with a whl sample
  int start_res_value; // starting time of a whl sample in res value
  int end_res_value; // end time of a whl sample
  int res_value_to_add; // number of res inside the interval
  int interval_index=0;
  double time_to_add;
  int bin_x; // where to add the time in the map
  int bin_y;
  // set the map to 0
  for(int x = 0; x < x_bins; x++)
    for(int y = 0; y < y_bins; y++)
      map[(x*y_bins) + y] = 0;

  // loop with the whl data
  for(int i = 0; i < whl_lines; i++)
    {
      if (x_whl[i] != -1.0)
	{
	  // get the begining of the period cover by this whl sample, in res value
	  res_value_for_whl_sample=(res_samples_per_whl_sample*i)+res_samples_per_whl_sample;
	  start_res_value=res_value_for_whl_sample-(res_samples_per_whl_sample/2); // start of whl sample
	  end_res_value=res_value_for_whl_sample+(res_samples_per_whl_sample/2); // end of whl sample
	  // start with the assumption that time to add is  = 0; as if outside the valid intervals
	  res_value_to_add=0;
	  
	  // then add res_time if there are intervals between start_res_value and end_res_value
	  while(start_interval[interval_index]<end_res_value && interval_index<interval_lines)
	    {
	      if(end_interval[interval_index]>start_res_value) // need to add some time
		{
		  // add the total interval time, and then remove what is not overlapping the whl sample
		  res_value_to_add=end_interval[interval_index]-start_interval[interval_index]; 
		  
		  if (start_interval[interval_index]<start_res_value)
		    {
		      res_value_to_add=res_value_to_add-(start_res_value-start_interval[interval_index]);
		    }
		  if (end_interval[interval_index]>end_res_value)
		    {
		      res_value_to_add=res_value_to_add-(end_interval[interval_index]-end_res_value);
		    }
		}
	      interval_index++; // go to the next interval
	    }
	  interval_index--;
	  // calculate the time in ms, and add it to the right bin of the occ map
	  time_to_add=ms_per_sample*((double)res_value_to_add/res_samples_per_whl_sample);
	  bin_x = (int) x_whl[i]/pixels_per_bin_x;
	  bin_y = (int) y_whl[i]/pixels_per_bin_y;

	  // check if it falls in the map
	  if ((bin_x<x_bins)&&(bin_y<y_bins))
	    {
	      map[(bin_x * y_bins) + bin_y]=
		map[(bin_x * y_bins) + bin_y] +
		time_to_add;
	    }
	}
    }
  
  // set the bins that have less then min_occ ms
  for(int x = 0; x < x_bins; x++)
    {
      for(int y = 0; y < y_bins; y++)
	{
	  if (map[(x*y_bins) + y] <= min_occ)
	    {
	      map[(x*y_bins) + y] = -1.0;
	    }
	}
    }
  // remove the bins that are isolated in space
  int moving_x;
  int moving_y;
  int factor=1;
  int smoothing_line = (factor*2)+1;
  int moving_index;
  int adjacent_pos;
  int index;
  for(int x = 0; x < x_bins; x++)
    {
      for(int y = 0; y < y_bins; y++)
	{
	  index = (x*y_bins) + y;
	  if (map[index] != -1.0)
	    {
	      adjacent_pos=0;
 	      // move to top left and then add valid bin to data
	      moving_x = x - factor;
	      for(int xx = 0; xx < smoothing_line; xx++)
		{
		  moving_y = y - factor;
		  for(int yy = 0; yy < smoothing_line; yy++)
		    {
		      if((moving_x >= 0) && 
			 (moving_x < x_bins) &&
			 (moving_y >= 0) &&
			 (moving_y < y_bins))
			{
 			  moving_index = (moving_x * y_bins) + moving_y;
			  
			  if ( map[moving_index] > 0)
			    {adjacent_pos++;}
			}
		      moving_y++;
		    }
		  moving_x++;
		}
	      if (adjacent_pos < 2)
		{
		  map[index] = -1.0;
		}
	    }
	}
    }
}
SEXP firing_rate_map_2d_cwrap(SEXP num_bins_x_r,
			      SEXP num_bins_y_r, 
			      SEXP pixels_per_bin_x_r,
			      SEXP pixels_per_bin_y_r,
			      SEXP x_spike_r, 
			      SEXP y_spike_r,
			      SEXP clu_r,
			      SEXP res_lines_r,
			      SEXP cells_r,
			      SEXP num_cells_r,
			      SEXP occ_map_r,
			      SEXP smooth_map_sd_r)
{

  int num_bins_x=INTEGER_VALUE(num_bins_x_r);
  int num_bins_y=INTEGER_VALUE(num_bins_y_r);
  int num_cells=INTEGER_VALUE(num_cells_r);
  int total_size = num_bins_x*num_bins_y*num_cells;
  int* cells = INTEGER_POINTER(cells_r);
  int target_cell;
  SEXP out = PROTECT(allocVector(REALSXP,total_size));
  double* maps = (double*)malloc(total_size*sizeof(double)); 
  double* one_map;

  for(int i = 0; i < num_cells;i++)
    {
      target_cell = cells[i];
      // use a pointer to place give the address for each cell
      one_map = maps + (i*num_bins_x*num_bins_y);
      create_place_field(num_bins_x,
      			 num_bins_y,
      			 REAL(pixels_per_bin_x_r)[0],
      			 REAL(pixels_per_bin_y_r)[0],
      			 REAL(x_spike_r),
      			 REAL(y_spike_r),
      			 INTEGER_POINTER(clu_r),
      			 INTEGER_VALUE(res_lines_r),
      			 target_cell,
      			 REAL(occ_map_r),
      			 one_map);
      /// smooth the place field
      smooth_double_gaussian_2d(one_map,num_bins_x, num_bins_y, REAL(smooth_map_sd_r)[0],-1.0);
    }
  
  double* ans;
  ans=REAL(out);
  for(int i = 0; i < total_size;i++)
    ans[i]=maps[i];

  free(maps);
  UNPROTECT(1);
  return(out);
}

void create_place_field( int x_bins,
			 int y_bins, 
			 double pixels_per_bin_x,
			 double pixels_per_bin_y,
			 double *x_spike, 
			 double *y_spike,
			 int *clu,
			 int res_lines,
			 int target_cell,
			 double *occupancy_map,
			 double *place_field)
{
  /*******************************************************
    function to calculate the firing rate map of one cell
    needs an occupancy map and the x and y position
    of all the spikes
  *********************************************************/
  int total_bins = x_bins * y_bins;
  int bin_x;
  int bin_y;
  // set the place_field array to 0
  for (int i = 0; i < total_bins ; i++)
    {
      place_field[i] = 0;
    }
  
  // add the number of spikes in the place_field array
  for(int i = 0; i < res_lines; i++)
    {
      if (clu[i]==target_cell)
	{
	  if (x_spike[i] != -1.0)
	    {
	      bin_x = (int) x_spike[i]/pixels_per_bin_x;
	      bin_y = (int) y_spike[i]/pixels_per_bin_y;
	      // check that this is in the map
	      if ((bin_x<x_bins)&&(bin_y<y_bins))
		{
		  place_field[((bin_x * y_bins) + bin_y)]++;
		}
	    }	  
	}
    }
  ////////////// get the firing rate for every bin////////////////
  for (int i = 0; i < total_bins; i++)
    {
      if (occupancy_map[i] == -1.0)
	{
	  place_field[i] = -1.0;
	}

      if (occupancy_map[i] != -1.0)
	{
	  if (place_field[i] != 0) // 
	    {
	      place_field[i] = (double)place_field[i]/((double)occupancy_map[i]/1000);
	    }
	}
    }
}
SEXP information_score_cwrap(SEXP cells_r,
			     SEXP cell_lines_r,
			     SEXP all_rate_maps_r,
			     SEXP occupancy_map_r,
			     SEXP map_size_r)
{
  int cell_lines = INTEGER_VALUE(cell_lines_r);
  int map_size = INTEGER_VALUE(map_size_r);
  int* cells = INTEGER_POINTER(cells_r);
  double* maps = REAL(all_rate_maps_r);
  double* occ_map = REAL(occupancy_map_r);
  double* one_map;
  SEXP out = PROTECT(allocVector(REALSXP,cell_lines));
  double* info = REAL(out);
  

  for(int i = 0; i < cell_lines; i++){
    one_map=maps+(i*map_size);
    info[i]=information_score(one_map,occ_map,map_size);
  }

  UNPROTECT(1);
  return(out);
}



double information_score(double* fr_map,
			 double* occ_map,
			 int map_size)
{
  
  /////////calculation for information score ////////////////
  /// see Skaggs et al. 1996, hippocampus
  double information_score=0;
  double total_occupancy=0;
  double* occupancy_probability;
  double info_j;
  double info_h;
  double global_fr;
  occupancy_probability = (double*)malloc(map_size*sizeof(double));
  // sum occupancy map to create occupancy probability
  for (int i = 0; i < map_size; i++)
    {
      if (occ_map[i]!=-1.0)
	{
	  total_occupancy += occ_map[i];
	} 
    }
  for (int i = 0; i < map_size; i++)
    {
      if (occ_map[i] != -1.0)
	{
	  occupancy_probability[i] = occ_map[i]/ total_occupancy;
  	} 
    }
  // calculate mean firing rate 
  global_fr=0;
  int valid = 0;
  for (int i = 0; i < map_size; i++)
    {
      if (occ_map[i] != -1.0)
	{
	  valid++;
	  global_fr+=fr_map[i]*occupancy_probability[i];
	}
    }

  // information score formula
  for (int i = 0; i < map_size; i++)
    {
      if(occ_map[i] != -1.0 && fr_map[i] !=-1.0)
	{
	  if (fr_map[i] == 0)
	    {
	      info_j = 0; 
	      info_h = 0;
	    }
	  else
	    { 
	      info_j = occupancy_probability[i]*fr_map[i]/global_fr;
	      info_h = log2(fr_map[i]/global_fr);
	    }
	  information_score += info_j*info_h;
	}
    }
  free(occupancy_probability);
  return information_score;
}

SEXP sparsity_score_cwrap(SEXP cells_r,
			  SEXP cell_lines_r,
			  SEXP all_rate_maps_r,
			  SEXP occupancy_map_r,
			  SEXP map_size_r)
{
  int cell_lines = INTEGER_VALUE(cell_lines_r);
  int map_size = INTEGER_VALUE(map_size_r);
  int* cells = INTEGER_POINTER(cells_r);
  double* maps = REAL(all_rate_maps_r);
  double* occ_map = REAL(occupancy_map_r);
  double* one_map;
  SEXP out = PROTECT(allocVector(REALSXP,cell_lines));
  double* info = REAL(out);
  

  for(int i = 0; i < cell_lines; i++){
    one_map=maps+(i*map_size);
    info[i]=sparsity_score(one_map,occ_map,map_size);
  }

  UNPROTECT(1);
  return(out);
}



double sparsity_score(double* fr_map,
		      double* occ_map,
		      int map_size)
{
  /////////calculation for sparsity score ////////////////
  /// see Skaggs et al. 1996, hippocampus
  /// reversed so that high score indicates high sparsity
  
  double sparsity_score=0;
  double total_occupancy=0;
  double* occupancy_probability;
  double sparsity_den=0;
  double sparsity_num=0;
  double sparsity_den_sum=0;
  double sparsity_num_sum=0;
  occupancy_probability = (double*)malloc(map_size*sizeof(double));
  // get a map with occupancy probability
  for (int i = 0; i < map_size; i++)
    {
      if (occ_map[i] == -1.0)
	{
	  occupancy_probability[i] = -1.0;
	}
      else
	{
	  occupancy_probability[i] = occ_map[i];
	  total_occupancy = total_occupancy + occ_map[i];
	} 
    }
  for (int i = 0; i < map_size; i++)
    {
      if (occ_map[i] != -1.0)
	{
	  occupancy_probability[i] = occupancy_probability[i]/ total_occupancy;
  	} 
    }

  for (int i = 0; i < map_size; i++)
    {
      if(occupancy_probability[i] != -1.0 && fr_map[i] !=-1.0)
	{
	  if (fr_map[i] == 0)
	    {
	      sparsity_num=0;
	      sparsity_den=0;
	    }
	  else
	    { 
	      sparsity_num=occupancy_probability[i]*fr_map[i];
	      sparsity_den=occupancy_probability[i]*(fr_map[i]*fr_map[i]);
	    }
	  sparsity_num_sum=sparsity_num_sum+sparsity_num;
	  sparsity_den_sum=sparsity_den_sum+sparsity_den;
	}
      sparsity_score=(sparsity_num_sum*sparsity_num_sum)/sparsity_den_sum;
    }
  free(occupancy_probability);
  return 1.0-sparsity_score;
}


SEXP map_autocorrelation_cwrap(SEXP cells_r,
			       SEXP cell_lines_r,
			       SEXP maps_r, 
			       SEXP num_bins_x_r,
			       SEXP num_bins_y_r,
			       SEXP auto_num_bins_x_r,
			       SEXP auto_num_bins_y_r,
			       SEXP min_bins_for_autocorrelation_r)
{
  int* cells = INTEGER_POINTER(cells_r);
  int cell_lines = INTEGER_VALUE(cell_lines_r);
  double* maps = REAL(maps_r);
  double* one_map;
  int num_bins_x = INTEGER_VALUE(num_bins_x_r);
  int num_bins_y = INTEGER_VALUE(num_bins_y_r);
  int auto_num_bins_x = INTEGER_VALUE(auto_num_bins_x_r);
  int auto_num_bins_y = INTEGER_VALUE(auto_num_bins_y_r);
  int min_bins_for_autocorrelation = INTEGER_VALUE(min_bins_for_autocorrelation_r);

  int total_bins_auto= auto_num_bins_x*auto_num_bins_y;
  SEXP out = PROTECT(allocVector(REALSXP,total_bins_auto*cell_lines));
  double* all_auto = REAL(out);
  double* one_auto;

  
  for(int i = 0; i < cell_lines; i++){
    one_auto=all_auto+(i*total_bins_auto);
    one_map=maps+(i*num_bins_x*num_bins_y);
    map_autocorrelation(one_map, // pointer to one place field map
			one_auto, // pointer to one spatial autocorrelation map
			num_bins_x, // x size of the place field map (num bins)
			num_bins_y, // y size of the place field map
			auto_num_bins_x, // x size of the autocorrelation map
			auto_num_bins_y, // y size of the autocorreltion map
			min_bins_for_autocorrelation); // minimum of valid values to do the correlation
  }

  UNPROTECT(1);
  return (out);
}


void map_autocorrelation(double *one_place, // pointer to one place field map
			 double *one_auto, // pointer to one spatial autocorrelation map
			 int x_bins_place_map, // x size of the place field map (num bins)
			 int y_bins_place_map, // y size of the place field map
			 int x_bins_auto_map, // x size of the autocorrelation map
			 int y_bins_auto_map, // y size of the autocorreltion map
			 int min_for_correlation) // minimum of valid values to do the correlation
{
/*************************************************************
 funciton to do the spatial autocorrelation for a place firing
map. 
 one_place should have a size = x_bins_place_map*y_bins_place_map
 x_bins_auto_map should = (x_bins_place_map*2)+1
 y_bins_auto_map should = (y_bins_place_map*2)+1
 one_auto should have a size =  x_bins_auto_map * y_bins_auto_map
*************************************************************/
  int min_x_offset = 0 - x_bins_place_map;
  int max_x_offset = x_bins_place_map;
  int min_y_offset = 0 - y_bins_place_map;
  int max_y_offset = y_bins_place_map;
  int mid_x = x_bins_place_map; // mid x value in autocorrelation
  int mid_y = y_bins_place_map; // min y value in autocorrelation
  int auto_x;
  int auto_y;
  int index_1;
  int index_2;
  int index_auto;
  int total_bins_place_map = x_bins_place_map * y_bins_place_map;
  int total_bins_auto_map = x_bins_auto_map * y_bins_auto_map;
  int offset_x;
  int offset_y;
 
  int n;
  double r;

  double* value_1_correlation;
  double* value_2_correlation;
  value_1_correlation = (double*)malloc(total_bins_place_map*sizeof(double));
  value_2_correlation = (double*)malloc(total_bins_place_map*sizeof(double));
  // set the auto_place map to -2, this is the invalid value, correlation range from -1 to 1
  for (int i = 0; i < total_bins_auto_map ; i++)
    {
      one_auto[i] = -2;
    }
  // loop for all possible lags in the x axis
  for (int x_off = min_x_offset; x_off <= max_x_offset; x_off++)
    {
      // loop for all possible lags in the y axis
      for (int y_off = min_y_offset; y_off <= max_y_offset; y_off++ )
	{
	  // for all the possible lags, calculate the following values
	  n = 0;  // number of valid lags to do the correlation for this offset
	  r = 0;  // r value of the correlation
	  // loop for all bins in the place map
	  for(int x = 0; x < x_bins_place_map; x++)
	    {
	      for(int y = 0; y < y_bins_place_map; y++)
		{
		  offset_x = x+x_off;
		  offset_y = y+y_off;
		  if ((offset_x >=0 && offset_x < x_bins_place_map) && (offset_y >= 0 && offset_y < y_bins_place_map))
		    {
		      index_1 = (x*y_bins_place_map) + y; // that is the index for the current bin in the place firing rate map
		      index_2 = ((offset_x)*y_bins_place_map)+(offset_y); // that is the index in the offset bin relative to the current bin
		      if (one_place[index_1]!=-1.0 &&one_place[index_2]!=-1.0) // -1 is the invalid value in the place firing rate map, only take into account the data if not invalid value 
			{
			  // copy the value into 2 vectors for the correlation
			  value_1_correlation[n]=one_place[index_1];
			  value_2_correlation[n]=one_place[index_2];
			  n++; 
			}
		    }
		}   
	    }
	  // if enough valid data to calculate the r value, if not the value for this lag will stay -2 
	  if ( n > min_for_correlation)
	    {
	      // calculate a correlation
	      r = correlation(value_1_correlation,value_2_correlation,n,-1.0);
	      auto_x = mid_x + x_off;
	      auto_y = mid_y + y_off;
	      index_auto= (auto_x*y_bins_auto_map)+auto_y;
	      one_auto[index_auto]=r;
	    }
	}
    }
  free(value_1_correlation);
  free(value_2_correlation);
}
