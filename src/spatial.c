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
		    double px_per_cm,
		    int res_sampling_rate, // ex 20 000
		    int res_samples_per_whl_sample) // ex 512
{
  double samples_per_sec=(double)res_sampling_rate/(double)res_samples_per_whl_sample;
  double diff_x;
  double diff_y;
  double distance_1;
  double distance_2;
  double distance_cm;
  // the speed at index x is estimated from position x-1 to position x+1
  // this is done to avoid introducing a shift, for example speed at index x is estimated from x and x+1
  // speed for the first and last samples is set to -1.0
  // if one of the 3 positions has -1.0, speed == -1.0
  
  speed[0]=-1.0;
  speed[whl_lines-1]=-1.0;
  
  for (int i = 1; i < whl_lines-1; i++)
  {
    if (x_whl[i-1] == -1.0 || y_whl[i-1] == -1.0 ||
        x_whl[i] == -1.0 || y_whl[i] == -1.0|| 
        x_whl[i+1]==-1.0 || y_whl[i+1]==-1.0)
    {
      speed[i]=-1.0;
    }
    else
    {
      // get the distance between previous-current, and current-after
      diff_x = sqrt(pow(x_whl[i-1] - x_whl[i],2)); 
      diff_y = sqrt(pow(y_whl[i-1] - y_whl[i],2));
      distance_1= sqrt(pow(diff_x,2)+pow(diff_y,2));
      
      diff_x = sqrt(pow(x_whl[i] - x_whl[i+1],2)); 
      diff_y = sqrt(pow(y_whl[i] - y_whl[i+1],2));
      distance_2= sqrt(pow(diff_x,2)+pow(diff_y,2));
      
      distance_cm= ((distance_1+distance_2)/2)/px_per_cm;
      
      speed[i] = distance_cm*samples_per_sec;
    }
  }
}


SEXP speed_from_whl_cwrap(SEXP x_whl_r,
			  SEXP y_whl_r,
			  SEXP whl_lines_r,
			  SEXP px_per_cm_r,
			  SEXP res_sampling_rate_r, 
			  SEXP res_samples_per_whl_sample_r)

{
  int whl_lines=INTEGER_VALUE(whl_lines_r);
  SEXP out = PROTECT(allocVector(REALSXP, whl_lines));
  double* speed = REAL(out);
  // call c function
  speed_from_whl(REAL(x_whl_r),
  		 REAL(y_whl_r),
  		 speed,
  		 whl_lines,
  		 REAL(px_per_cm_r)[0],
  		 INTEGER_VALUE(res_sampling_rate_r),
  		 INTEGER_VALUE(res_samples_per_whl_sample_r));

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
	      Rprintf("whl file is too small for some res values:%d\n",res[i]);
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
      Rprintf("speed_intervals_cwrap: number of intervals < 1");
      return(R_NilValue);
    }
  int* start = (int*)malloc(interval_lines*sizeof(int));
  int* end = (int*)malloc(interval_lines*sizeof(int));
  SEXP out = PROTECT(allocMatrix(REALSXP,interval_lines,2));
  speed_intervals(REAL(speed_r), 
		  whl_lines, 
		  INTEGER_VALUE(res_samples_per_whl_sample_r),
		  REAL(min_speed_r)[0],
		  REAL(max_speed_r)[0],
		  start,
		  end);
  double* rans;
  rans = REAL(out);
  for(int i = 0; i < interval_lines; i++) {
    rans[i] = start[i];
    rans[i+interval_lines]=end[i];
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
int direction_intervals_count(int* direction, int whl_lines, int res_samples_per_whl_sample, int target_direction)
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
    if (direction[i]==-1)
    {
      minus_one++;
      if (was_in==1 && minus_one > max_minus_one_before_out)
      { // to exit of interval when we lose track of animal
        end_int=(i*res_samples_per_whl_sample)+res_samples_per_whl_sample;
        was_in=0;
        jj=0;
        num_intervals++;
      }
    }
    if (direction[i]!=-1)
    {
      minus_one=0;
      if (was_in==1) // look for exit time and print
      {
        if (direction[i]!=target_direction)
        {
          end_int=(i*res_samples_per_whl_sample)+res_samples_per_whl_sample;
          jj=0;
          num_intervals++;
        }
      }
      if (was_in==0) // look for entry time in the speed interval
      {
        if(direction[i]==target_direction)
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
      }
    }
  }
  return(num_intervals);
}
void direction_intervals(int* direction, int whl_lines, int res_samples_per_whl_sample,int target_direction,int* start,int* end)
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
    if (direction[i]==-1)
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
      }
    }
    if (direction[i]!=-1)
    {
      minus_one=0;
      if (was_in==1) // look for exit time and print
      {
        if (direction[i]!=target_direction)
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
        if(direction[i]==target_direction)
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
      }
    }
  }
  return;
}

SEXP direction_intervals_cwrap(SEXP direction_r, SEXP whl_lines_r, SEXP res_samples_per_whl_sample_r, 
                               SEXP target_direction_r)
{
  int whl_lines=INTEGER_VALUE(whl_lines_r);
  int interval_lines;
  interval_lines=direction_intervals_count(INTEGER_POINTER(direction_r), 
                                            whl_lines,
                                            INTEGER_VALUE(res_samples_per_whl_sample_r),
                                            INTEGER_VALUE(target_direction_r));
  if(interval_lines<1)
  {
    Rprintf("direction_intervals_cwrap: number of intervals < 1");
    return(R_NilValue);
  }
  int* start = (int*)malloc(interval_lines*sizeof(int));
  int* end = (int*)malloc(interval_lines*sizeof(int));
  SEXP out = PROTECT(allocMatrix(REALSXP,interval_lines,2));
  direction_intervals(INTEGER_POINTER(direction_r), 
                      whl_lines, 
                      INTEGER_VALUE(res_samples_per_whl_sample_r),
                      INTEGER_VALUE(target_direction_r),
                      start,
                      end);
  double* rans;
  rans = REAL(out);
  for(int i = 0; i < interval_lines; i++) {
    rans[i] = start[i];
    rans[i+interval_lines]=end[i];
  }
  free(start);
  free(end);
  UNPROTECT(1);
  return(out);
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
                          end_interval_index,
                          1);
  
 // for(int i = 0; i < interval_lines;i++)
 //   Rprintf("int %d, start: %d, end: %d\n",i,start_interval_index[i],end_interval_index[i]);
  // set all the spike position to -1
  set_array_to_value_double(x_spike,res_lines,-1.0);
  set_array_to_value_double(y_spike,res_lines,-1.0);
  
  //loop for each interval and set the position when valid//
  int whl_index;
  int residual;
  for (int k = 0; k < interval_lines; k++)
  {
    for (int i = start_interval_index[k]; i <= end_interval_index[k]; i++)
    {
      if (res[i] != -1.0)
      {
        whl_index = (res[i] / res_samples_per_whl_sample) - 1;
        if (whl_index >=  whl_lines || (whl_index < 0 && whl_index != -1))
        {
          Rprintf("whl file is too small for some res value\n");
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
  SEXP out = PROTECT(allocMatrix(REALSXP,x_bins,y_bins)); // nrow ncol, to return, y and x are swapped 
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
      rans[y*x_bins+x] = map[x*y_bins+y]; // the way things are ordered in C differs from R
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
  int one_map_size=num_bins_x*num_bins_y;
  int* cells = INTEGER_POINTER(cells_r);
  int target_cell;
  double* occ=REAL(occ_map_r);
  SEXP out = PROTECT(allocVector(REALSXP,total_size));
  double* maps = (double*)malloc(total_size*sizeof(double)); 
  double* t_occ_map = (double*)malloc(one_map_size*sizeof(double));
  double* one_map;

  
  // transpose the occ maps from R to C format
  for(int x =0; x < num_bins_x;x++)
    for(int y =0; y < num_bins_y;y++)
      t_occ_map[x*num_bins_y+y]=occ[x+num_bins_x*y];
  
  
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
      			 t_occ_map,
      			 one_map);
      /// smooth the place field
      smooth_double_gaussian_2d(one_map,num_bins_x, num_bins_y, REAL(smooth_map_sd_r)[0],-1.0);
    }
  
  double* ans;
  double* one_ans;
  ans=REAL(out);
  
  for(int i = 0; i < total_size;i++){
    ans[i]=maps[i];
  }
  //transpose the maps for R
  for(int i = 0; i < num_cells;i++){
    one_map = maps + (i*one_map_size);
    one_ans = ans + (i*one_map_size);
    for(int x =0; x < num_bins_x;x++)
      for(int y =0; y < num_bins_y;y++)
        one_ans[x+num_bins_x*y]=one_map[x*num_bins_y+y];
  }

  free(maps);
  free(t_occ_map);
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
SEXP information_score_cwrap(SEXP cell_lines_r,
			     SEXP all_rate_maps_r,
			     SEXP occupancy_map_r,
			     SEXP map_size_r)
{
  int cell_lines = INTEGER_VALUE(cell_lines_r);
  int map_size = INTEGER_VALUE(map_size_r);
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

SEXP sparsity_score_cwrap(SEXP cell_lines_r,
			  SEXP all_rate_maps_r,
			  SEXP occupancy_map_r,
			  SEXP map_size_r)
{
  int cell_lines = INTEGER_VALUE(cell_lines_r);
  int map_size = INTEGER_VALUE(map_size_r);
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

SEXP map_autocorrelation_cwrap(SEXP cell_lines_r,
			       SEXP maps_r, 
			       SEXP num_bins_x_r,
			       SEXP num_bins_y_r,
			       SEXP auto_num_bins_x_r,
			       SEXP auto_num_bins_y_r,
			       SEXP min_bins_for_autocorrelation_r)
{
  int cell_lines = INTEGER_VALUE(cell_lines_r);
  double* maps = REAL(maps_r);
  double* one_map;
  double* one_tmap;
  int num_bins_x = INTEGER_VALUE(num_bins_x_r);
  int num_bins_y = INTEGER_VALUE(num_bins_y_r);
  int one_map_size=num_bins_x*num_bins_y;
  int auto_num_bins_x = INTEGER_VALUE(auto_num_bins_x_r);
  int auto_num_bins_y = INTEGER_VALUE(auto_num_bins_y_r);
  int min_bins_for_autocorrelation = INTEGER_VALUE(min_bins_for_autocorrelation_r);

  // transpose the maps
  double* tmaps = (double*) malloc(one_map_size*cell_lines*sizeof(double));
    //transpose the maps for R
  for(int i = 0; i < cell_lines;i++){
    one_map = maps + (i*one_map_size);
    one_tmap = tmaps + (i*one_map_size);
      for(int x =0; x < num_bins_x;x++)
        for(int y =0; y < num_bins_y;y++)
          one_tmap[x*num_bins_y+y]=one_map[x+num_bins_x*y];
    }
  
  int total_bins_auto= auto_num_bins_x*auto_num_bins_y;
  SEXP out = PROTECT(allocVector(REALSXP,total_bins_auto*cell_lines));
  double* all_auto = (double*) malloc(total_bins_auto*cell_lines*sizeof(double));
  double* one_auto;
  double* all_tauto=REAL(out);
  double* one_tauto;

  for(int i = 0; i < cell_lines; i++){
    one_auto=all_auto+(i*total_bins_auto);
    one_map=tmaps+(i*num_bins_x*num_bins_y);
    map_autocorrelation(one_map, // pointer to one place field map
			one_auto, // pointer to one spatial autocorrelation map
			num_bins_x, // x size of the place field map (num bins)
			num_bins_y, // y size of the place field map
			auto_num_bins_x, // x size of the autocorrelation map
			auto_num_bins_y, // y size of the autocorreltion map
			min_bins_for_autocorrelation); // minimum of valid values to do the correlation
  }

  //for(int i = 0; i < total_bins_auto*cell_lines;i++)
  //    all_tauto[i]=all_auto[i];
  
  
  // transpose the results
  for(int i = 0; i < cell_lines;i++){
    one_auto = all_auto + (i*total_bins_auto);
    one_tauto = all_tauto + (i*total_bins_auto);
    for(int x =0; x < auto_num_bins_x;x++)
      for(int y =0; y < auto_num_bins_y;y++)
        one_tauto[x+auto_num_bins_x*y]=one_auto[x*auto_num_bins_y+y];
  }

  
  free(tmaps);
  free(all_auto);
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


void get_x_and_y_bin_from_index(int x_bins,
				int y_bins,
				int index,
				int* x_bin,
				int* y_bin)
{
  //index=((bin_x * y_bins) + bin_y)
  *x_bin=index/y_bins;
  *y_bin=index%y_bins;
}
int get_index_from_x_and_y_bin(int x_bins,
			       int y_bins,
			       int x_bin,
			       int y_bin)
{
  return ((x_bin * y_bins) + y_bin);
}

void detect_one_field_with_field( double* map,
		      int x_bins,
		      int y_bins,
		      int min_num_bins_fields,
		      double threshold,
		      double* mean_x_field,
		      double* mean_y_field,
		      double* max_radius_field,
		      int* num_bins_field,
		      double invalid,
		      double* field)
{
  // function to detect on field starting from the peak value
  // the pixels for which the function tries to cumulate to form a field are set to invalid
  // so that the function do not always work with the same pixels if called many times
  // should be a recursive function but have no time to re-write
  int index;
  int* x_positive_pixels; // x value of the positive pixels
  int* y_positive_pixels; // y value of the positive pixels
  int num_positive_pixels=0; 
  int num_bins= x_bins*y_bins; // in the map
  int max_index; // index of the peak in the map
  double max; // peak value
  int x_max; // x value of the peak
  int y_max; // y value of the peak
  int margin; // when distance around the peak to add pixels
  int x_check;
  int y_check;
  int x_adj; // to find positive adjacent pixels
  int y_adj; 
  int added; // flag to see if a pixel was added for a given margin, true:1 false:0
  int* map_positive; // flag for pixels that have been added to field
  double possible_radius;
  map_positive = (int*)malloc(num_bins*sizeof(int));
  x_positive_pixels = (int*)malloc(num_bins*sizeof(int));
  y_positive_pixels = (int*)malloc(num_bins*sizeof(int));
  for (int i = 0 ; i < num_bins; i++)
    {
      map_positive[i]=0;
    }
  // set the field array to invalid
  set_array_to_value_double(field,num_bins,invalid);
  // find the max value in the map
  max=find_max_double_index(num_bins,map,&max_index);
  if(max<threshold)
    {  
      *mean_x_field=-1;
      *mean_y_field=-1;
      *max_radius_field=-1;
    }
  else
    { // loop from the max pixel and try to add as many positive pixels as possible
      num_positive_pixels=0;
      get_x_and_y_bin_from_index(x_bins,
				 y_bins,
				 max_index,
				 &x_max,
				 &y_max);
      x_positive_pixels[num_positive_pixels]= x_max;
      y_positive_pixels[num_positive_pixels]= y_max;
      field[max_index]=map[max_index];
      map[max_index]= invalid;
      map_positive[max_index]=1;
      margin=0;
      added=1;
      // to be added a pixel needs to  1) be above threshold and
      // 2) contact with a positive pixel
      while (added==1)
	{
	  added=0;
	  margin++;
	  //
	  // do the left vertical band
	  //
	  x_check=x_max-margin;
	  
	  if(x_check>=0)
	    {
	      for (y_check=y_max-margin;y_check <= y_max+margin; y_check++)
		{
		  if (y_check >=0 && y_check < y_bins) // if y within legal range
		    {
		      if (map[get_index_from_x_and_y_bin(x_bins,y_bins,x_check,y_check)] > threshold)
			{
			  // we found a pixel in the margin that is above threshold, now look 
			  // if adjacent pixels were set true in map_positive 
			  for (x_adj = x_check -1 ; x_adj <= x_check + 1 ;x_adj++)
			    {
			      if (x_adj>=0 && x_adj < x_bins)
				for (y_adj = y_check -1; y_adj <= y_check+1; y_adj++)
				  {
				    if (y_adj>=0 && y_adj < y_bins)
				      {
					if (map_positive[get_index_from_x_and_y_bin(x_bins,y_bins,x_adj,y_adj)] == 1)
					  { // we found a pixel adjacent the pixel of interest that is positive
					    // set the current pixel positive
					    added=1;
					    map_positive[get_index_from_x_and_y_bin(x_bins,y_bins,x_check,y_check)]=1;
					    num_positive_pixels++;
					    x_positive_pixels[num_positive_pixels]=x_check;
					    y_positive_pixels[num_positive_pixels]=y_check;
					    index=get_index_from_x_and_y_bin(x_bins,
									     y_bins,
									     x_positive_pixels[num_positive_pixels],
									     y_positive_pixels[num_positive_pixels]);
					    field[index]=map[index];
					    map[index]= invalid;
					    // two for loops because the pixel was already added ////
					    y_adj=y_check+1;
					    x_adj=x_check+1;
					  }
				      }
				  }
			    }
			}
		      
		    }
		}
	    }
	  //
	  // do the right vertical band
	  //
	  
	  x_check=x_max+margin;
	  if(x_check<x_bins)
	    {
	      for (y_check=y_max-margin;y_check <= y_max+margin; y_check++)
		{
		  if (y_check >=0 && y_check < y_bins)
		    {
		      if (map[get_index_from_x_and_y_bin(x_bins,y_bins,x_check,y_check)] > threshold)
			{
			  // we found a pixel in the margin that is above threshold, now look 
			  // if adjacent pixels were set true in map_positive 
			  for (x_adj = x_check -1 ; x_adj <= x_check + 1 ;x_adj++)
			    {
			      if (x_adj>=0 && x_adj < x_bins)
				for (y_adj = y_check -1; y_adj <= y_check+1; y_adj++)
				  {
				    if (y_adj>=0 && y_adj < y_bins)
				      {
					if (map_positive[get_index_from_x_and_y_bin(x_bins,y_bins,x_adj,y_adj)] == 1)
					  { // we found a pixel adjacent the pixel of interest that is positive
					    // set the current pixel positive
					    added=1;
					    map_positive[get_index_from_x_and_y_bin(x_bins,y_bins,x_check,y_check)]=1;
					    num_positive_pixels++;
					    x_positive_pixels[num_positive_pixels]=x_check;
					    y_positive_pixels[num_positive_pixels]=y_check;
					    index=get_index_from_x_and_y_bin(x_bins,
									     y_bins,
									     x_positive_pixels[num_positive_pixels],
									     y_positive_pixels[num_positive_pixels]);
					    field[index]=map[index];
					    map[index]= invalid;
					    // should find a way to exit the two for loop to save time ////
					    y_adj=y_check+1;
					    x_adj=x_check+1;
					  }
				      }
				  }
			    }
			}
		      
		    }
		}
	    }
	  //
	  // do the top horizontal band
	  //
	  y_check=y_max+margin;
	  if(y_check<y_bins)
	    {
	      for (x_check=x_max-margin+1;x_check <= x_max+margin-1; x_check++) 
		{
		  if (x_check >=0 && x_check < x_bins)
		    {
		      if (map[get_index_from_x_and_y_bin(x_bins,y_bins,x_check,y_check)] > threshold)
			{
			  // we found a pixel in the margin that is above threshold, now look 
			  // if adjacent pixels were set true in map_positive 
			  for (x_adj = x_check -1 ; x_adj <= x_check + 1 ;x_adj++)
			    {
			      if (x_adj>=0 && x_adj < x_bins)
				for (y_adj = y_check -1; y_adj <= y_check+1; y_adj++)
				  {
				    if (y_adj>=0 && y_adj < y_bins)
				      {
					if (map_positive[get_index_from_x_and_y_bin(x_bins,y_bins,x_adj,y_adj)] == 1)
					  { // we found a pixel adjacent the pixel of interest that is positive
					    // set the current pixel positive
					    added=1;
					    map_positive[get_index_from_x_and_y_bin(x_bins,y_bins,x_check,y_check)]=1;
					    num_positive_pixels++;
					    x_positive_pixels[num_positive_pixels]=x_check;
					    y_positive_pixels[num_positive_pixels]=y_check;
					    index=get_index_from_x_and_y_bin(x_bins,
									     y_bins,
									     x_positive_pixels[num_positive_pixels],
									     y_positive_pixels[num_positive_pixels]);
					    field[index]=map[index];
					    map[index]= invalid;
					    // should find a way to exit the two for loop to save time ////
					    y_adj=y_check+1;
					    x_adj=x_check+1;
					  }
				      }
				  }
			    }
			}
		      
		    }
		}
	    }

	  //
	  // do the bottom horizontal band
	  //
	  y_check=y_max-margin;
	  if(y_check>=0)
	    {
	      for (x_check=x_max-margin+1;x_check <= x_max+margin-1; x_check++) 
		{
		  if (x_check >=0 && x_check < x_bins)
		    {
		      if (map[get_index_from_x_and_y_bin(x_bins,y_bins,x_check,y_check)] > threshold)
			{
			  // we found a pixel in the margin that is above threshold, now look 
			  // if adjacent pixels were set true in map_positive 
			  for (x_adj = x_check -1 ; x_adj <= x_check + 1 ;x_adj++)
			    {
			      if (x_adj>=0 && x_adj < x_bins)
				for (y_adj = y_check -1; y_adj <= y_check+1; y_adj++)
				  {
				    if (y_adj>=0 && y_adj < y_bins)
				      {
					if (map_positive[get_index_from_x_and_y_bin(x_bins,y_bins,x_adj,y_adj)] == 1)
					  { // we found a pixel adjacent the pixel of interest that is positive
					    // set the current pixel positive
					    added=1;
					    map_positive[get_index_from_x_and_y_bin(x_bins,y_bins,x_check,y_check)]=1;
					    num_positive_pixels++;
					    x_positive_pixels[num_positive_pixels]=x_check;
					    y_positive_pixels[num_positive_pixels]=y_check;
					    index=get_index_from_x_and_y_bin(x_bins,
									     y_bins,
									     x_positive_pixels[num_positive_pixels],
									     y_positive_pixels[num_positive_pixels]);
					    field[index]=map[index];
					    map[index]= invalid;
					    // should find a way to exit the two for loop to save time ////
					    y_adj=y_check+1;
					    x_adj=x_check+1;
					  }
				      }
				  }
			    }
			}
		      
		    }
		}
	    }
	}
      // if enough pixels were detected to form a field, get the mean x and y 
      // and set the value to invalid in map
      if (num_positive_pixels >= min_num_bins_fields)
	{
	  *num_bins_field=num_positive_pixels+1;// because it is not the index now
	  *mean_x_field=0;
	  *mean_y_field=0;
	  *max_radius_field=0;
	  for (int i = 0; i < num_positive_pixels; i++)
	    {
	      *mean_x_field=(*mean_x_field)+x_positive_pixels[i];
	      *mean_y_field=(*mean_y_field)+y_positive_pixels[i];
	    }
	  *mean_x_field=(*mean_x_field)/num_positive_pixels;
	  *mean_y_field=(*mean_y_field)/num_positive_pixels;
	  for (int i = 0; i < num_positive_pixels; i++)
	    {
	      possible_radius=distance(*mean_x_field,*mean_y_field,x_positive_pixels[i],y_positive_pixels[i]);
	      if (possible_radius > (*max_radius_field))
		{
		  *max_radius_field=possible_radius;
		}
	    }
	}
      else
	{

	  *mean_x_field=-1;
	  *mean_y_field=-1;
	  *max_radius_field=-1;
	  *num_bins_field=-1;
	}
    }
  free(y_positive_pixels);
  free(x_positive_pixels);
  free(map_positive);
  return;
}


void detect_one_field( double* map,
		      int x_bins,
		      int y_bins,
		      int min_num_bins_fields,
		      double threshold,
		      double* mean_x_field,
		      double* mean_y_field,
		      double* max_radius_field,
		      int* num_bins_field,
		      double invalid)
{
  // function to detect on field starting from the peak value
  // the pixels for which the function tries to cumulate to form a field are set to invalid
  // so that the function do not always work with the same pixels if called many times
  // should be a recursive function but have no time to re-write
  int index;
  int* x_positive_pixels; // x value of the positive pixels
  int* y_positive_pixels; // y value of the positive pixels
  int num_positive_pixels=0; 
  int num_bins= x_bins*y_bins; // in the map
  int max_index; // index of the peak in the map
  double max; // peak value
  int x_max; // x value of the peak
  int y_max; // y value of the peak
  int margin; // when distance around the peak to add pixels
  int x_check;
  int y_check;
  int x_adj; // to find positive adjacent pixels
  int y_adj; 
  int added; // flag to see if a pixel was added for a given margin, true:1 false:0
  int* map_positive; // flag for pixels that have been added to field
  double possible_radius;
  map_positive = (int*)malloc(num_bins*sizeof(int));
  x_positive_pixels = (int*)malloc(num_bins*sizeof(int));
  y_positive_pixels = (int*)malloc(num_bins*sizeof(int));
  for (int i = 0 ; i < num_bins; i++)
    {
      map_positive[i]=0;
    }

  // find the max value in the map
  max=find_max_double_index(num_bins,map,&max_index);
  if(max<threshold)
    {  
      *mean_x_field=-1;
      *mean_y_field=-1;
      *max_radius_field=-1;
    }
  else
    { // loop from the max pixel and try to add as many positive pixels as possible
      num_positive_pixels=0;
      get_x_and_y_bin_from_index(x_bins,
				 y_bins,
				 max_index,
				 &x_max,
				 &y_max);
      x_positive_pixels[num_positive_pixels]= x_max;
      y_positive_pixels[num_positive_pixels]= y_max;
      map[max_index]= invalid;
      map_positive[max_index]=1;
      margin=0;
      added=1;
      // to be added a pixel needs to  1) be above threshold and
      // 2) contact with a positive pixel
      while (added==1)
	{
	  added=0;
	  margin++;
	  //
	  // do the left vertical band
	  //
	  x_check=x_max-margin;
	  
	  if(x_check>=0)
	    {
	      for (y_check=y_max-margin;y_check <= y_max+margin; y_check++)
		{
		  if (y_check >=0 && y_check < y_bins) // if y within legal range
		    {
		      if (map[get_index_from_x_and_y_bin(x_bins,y_bins,x_check,y_check)] > threshold)
			{
			  // we found a pixel in the margin that is above threshold, now look 
			  // if adjacent pixels were set true in map_positive 
			  for (x_adj = x_check -1 ; x_adj <= x_check + 1 ;x_adj++)
			    {
			      if (x_adj>=0 && x_adj < x_bins)
				for (y_adj = y_check -1; y_adj <= y_check+1; y_adj++)
				  {
				    if (y_adj>=0 && y_adj < y_bins)
				      {
					if (map_positive[get_index_from_x_and_y_bin(x_bins,y_bins,x_adj,y_adj)] == 1)
					  { // we found a pixel adjacent the pixel of interest that is positive
					    // set the current pixel positive
					    added=1;
					    map_positive[get_index_from_x_and_y_bin(x_bins,y_bins,x_check,y_check)]=1;
					    num_positive_pixels++;
					    x_positive_pixels[num_positive_pixels]=x_check;
					    y_positive_pixels[num_positive_pixels]=y_check;
					    index=get_index_from_x_and_y_bin(x_bins,
									     y_bins,
									     x_positive_pixels[num_positive_pixels],
									     y_positive_pixels[num_positive_pixels]);
					    map[index]= invalid;
					    // two for loops because the pixel was already added ////
					    y_adj=y_check+1;
					    x_adj=x_check+1;
					  }
				      }
				  }
			    }
			}
		      
		    }
		}
	    }
	  //
	  // do the right vertical band
	  //
	  
	  x_check=x_max+margin;
	  if(x_check<x_bins)
	    {
	      for (y_check=y_max-margin;y_check <= y_max+margin; y_check++)
		{
		  if (y_check >=0 && y_check < y_bins)
		    {
		      if (map[get_index_from_x_and_y_bin(x_bins,y_bins,x_check,y_check)] > threshold)
			{
			  // we found a pixel in the margin that is above threshold, now look 
			  // if adjacent pixels were set true in map_positive 
			  for (x_adj = x_check -1 ; x_adj <= x_check + 1 ;x_adj++)
			    {
			      if (x_adj>=0 && x_adj < x_bins)
				for (y_adj = y_check -1; y_adj <= y_check+1; y_adj++)
				  {
				    if (y_adj>=0 && y_adj < y_bins)
				      {
					if (map_positive[get_index_from_x_and_y_bin(x_bins,y_bins,x_adj,y_adj)] == 1)
					  { // we found a pixel adjacent the pixel of interest that is positive
					    // set the current pixel positive
					    added=1;
					    map_positive[get_index_from_x_and_y_bin(x_bins,y_bins,x_check,y_check)]=1;
					    num_positive_pixels++;
					    x_positive_pixels[num_positive_pixels]=x_check;
					    y_positive_pixels[num_positive_pixels]=y_check;
					    index=get_index_from_x_and_y_bin(x_bins,
									     y_bins,
									     x_positive_pixels[num_positive_pixels],
									     y_positive_pixels[num_positive_pixels]);
					    map[index]= invalid;
					    // should find a way to exit the two for loop to save time ////
					    y_adj=y_check+1;
					    x_adj=x_check+1;
					  }
				      }
				  }
			    }
			}
		      
		    }
		}
	    }
	  //
	  // do the top horizontal band
	  //
	  y_check=y_max+margin;
	  if(y_check<y_bins)
	    {
	      for (x_check=x_max-margin+1;x_check <= x_max+margin-1; x_check++) 
		{
		  if (x_check >=0 && x_check < x_bins)
		    {
		      if (map[get_index_from_x_and_y_bin(x_bins,y_bins,x_check,y_check)] > threshold)
			{
			  // we found a pixel in the margin that is above threshold, now look 
			  // if adjacent pixels were set true in map_positive 
			  for (x_adj = x_check -1 ; x_adj <= x_check + 1 ;x_adj++)
			    {
			      if (x_adj>=0 && x_adj < x_bins)
				for (y_adj = y_check -1; y_adj <= y_check+1; y_adj++)
				  {
				    if (y_adj>=0 && y_adj < y_bins)
				      {
					if (map_positive[get_index_from_x_and_y_bin(x_bins,y_bins,x_adj,y_adj)] == 1)
					  { // we found a pixel adjacent the pixel of interest that is positive
					    // set the current pixel positive
					    added=1;
					    map_positive[get_index_from_x_and_y_bin(x_bins,y_bins,x_check,y_check)]=1;
					    num_positive_pixels++;
					    x_positive_pixels[num_positive_pixels]=x_check;
					    y_positive_pixels[num_positive_pixels]=y_check;
					    index=get_index_from_x_and_y_bin(x_bins,
									     y_bins,
									     x_positive_pixels[num_positive_pixels],
									     y_positive_pixels[num_positive_pixels]);
					    map[index]= invalid;
					    // should find a way to exit the two for loop to save time ////
					    y_adj=y_check+1;
					    x_adj=x_check+1;
					  }
				      }
				  }
			    }
			}
		      
		    }
		}
	    }

	  //
	  // do the bottom horizontal band
	  //
	  y_check=y_max-margin;
	  if(y_check>=0)
	    {
	      for (x_check=x_max-margin+1;x_check <= x_max+margin-1; x_check++) 
		{
		  if (x_check >=0 && x_check < x_bins)
		    {
		      if (map[get_index_from_x_and_y_bin(x_bins,y_bins,x_check,y_check)] > threshold)
			{
			  // we found a pixel in the margin that is above threshold, now look 
			  // if adjacent pixels were set true in map_positive 
			  for (x_adj = x_check -1 ; x_adj <= x_check + 1 ;x_adj++)
			    {
			      if (x_adj>=0 && x_adj < x_bins)
				for (y_adj = y_check -1; y_adj <= y_check+1; y_adj++)
				  {
				    if (y_adj>=0 && y_adj < y_bins)
				      {
					if (map_positive[get_index_from_x_and_y_bin(x_bins,y_bins,x_adj,y_adj)] == 1)
					  { // we found a pixel adjacent the pixel of interest that is positive
					    // set the current pixel positive
					    added=1;
					    map_positive[get_index_from_x_and_y_bin(x_bins,y_bins,x_check,y_check)]=1;
					    num_positive_pixels++;
					    x_positive_pixels[num_positive_pixels]=x_check;
					    y_positive_pixels[num_positive_pixels]=y_check;
					    index=get_index_from_x_and_y_bin(x_bins,
									     y_bins,
									     x_positive_pixels[num_positive_pixels],
									     y_positive_pixels[num_positive_pixels]);
					    map[index]= invalid;
					    // should find a way to exit the two for loop to save time ////
					    y_adj=y_check+1;
					    x_adj=x_check+1;
					  }
				      }
				  }
			    }
			}
		      
		    }
		}
	    }
	}
      // if enough pixels were detected to form a field, get the mean x and y 
      // and set the value to invalid in map
      if (num_positive_pixels >= min_num_bins_fields)
	{
	  *num_bins_field=num_positive_pixels+1;// because it is not the index now
	  *mean_x_field=0;
	  *mean_y_field=0;
	  *max_radius_field=0;
	  for (int i = 0; i < num_positive_pixels; i++)
	    {
	      *mean_x_field=(*mean_x_field)+x_positive_pixels[i];
	      *mean_y_field=(*mean_y_field)+y_positive_pixels[i];
	    }
	  *mean_x_field=(*mean_x_field)/num_positive_pixels;
	  *mean_y_field=(*mean_y_field)/num_positive_pixels;
	  for (int i = 0; i < num_positive_pixels; i++)
	    {
	      possible_radius=distance(*mean_x_field,*mean_y_field,x_positive_pixels[i],y_positive_pixels[i]);
	      if (possible_radius > (*max_radius_field))
		{
		  *max_radius_field=possible_radius;
		}
	    }
	}
      else
	{

	  *mean_x_field=-1;
	  *mean_y_field=-1;
	  *max_radius_field=-1;
	  *num_bins_field=-1;
	}
    }
  free(y_positive_pixels);
  free(x_positive_pixels);
  free(map_positive);
  return;
}



void grid_closest_peaks_to_middle(double *map,
				  int x_bins,
				  int y_bins,
				  double pixels_per_bin,
				  int num_fields_to_detect,
				  int min_num_bins_fields,
				  float threshold,
				  double invalid,
				  double* x_field_res, // res for results
				  double* y_field_res,
				  double* radius_field_res,
				  double* distance_to_mid_res,
				  int* area_field_res,
				  int num_fields_in_res)
{
  /***********************************************************************
      detect the closest peaks to the middle of the map (include field
      at center)
  ************************************************************************/
  double* x_field;
  double* y_field;
  double* distance_to_mid;
  double* radius_field;
  int* area_field;
  int index=0;
  int num_valid_fields=0;
  int mid_x=x_bins/2;
  int mid_y=y_bins/2;

  x_field = (double*)malloc(num_fields_to_detect*sizeof(double));
  y_field = (double*)malloc(num_fields_to_detect*sizeof(double));
  distance_to_mid = (double*)malloc(num_fields_to_detect*sizeof(double));
  radius_field = (double*)malloc(num_fields_to_detect*sizeof(double));
  area_field = (int*)malloc(num_fields_to_detect*sizeof(int));
  // detect the fields
  for(int i = 0; i < num_fields_to_detect ; i++)
    { 
      detect_one_field(map,
		       x_bins,
		       y_bins,
		       min_num_bins_fields,
		       threshold,
		       & x_field[i],
		       & y_field[i],
		       & radius_field[i],
		       & area_field[i],
		       invalid);
      if(x_field[i]!=-1)
	{num_valid_fields++;}
    }
 
  set_array_to_value_double(x_field_res,num_fields_in_res,-1);
  set_array_to_value_double(y_field_res,num_fields_in_res,-1);
  set_array_to_value_double(radius_field_res,num_fields_in_res,-1);
  set_array_to_value_double(distance_to_mid_res,num_fields_in_res,-1);
  set_array_to_value_int(area_field_res,num_fields_in_res,-1);
  
  // get the distance between each field and center
  for(int i = 0; i < num_fields_to_detect ; i++)
    {
      if (x_field[i] != -1 && y_field[i] != -1)
	{
	  distance_to_mid[i]=distance(mid_x,mid_y,x_field[i],y_field[i]);
	}
      else
	{
	  distance_to_mid[i]=-1;
	}
    }
  // get the position of the fields closer to the center of the map
  for (int i = 0; i < num_fields_in_res && i < num_valid_fields; i++)
    {
      distance_to_mid_res[i]=999999999; // should be larger what is in the data
      for(int j = 0; j < num_fields_to_detect ; j++)
	{
	  if (distance_to_mid[j]!=-1)
	    {
	      if (distance_to_mid[j]<distance_to_mid_res[i])
		{
		  index=j;
		  // save the data for the results
		  distance_to_mid_res[i]=distance_to_mid[j];
		}
	    }
	}
      // save the results and set the smallest distance to -1
      x_field_res[i]=x_field[index]*pixels_per_bin;
      y_field_res[i]=y_field[index]*pixels_per_bin;
      radius_field_res[i]=radius_field[index]*pixels_per_bin;
      area_field_res[i]=area_field[index]*pixels_per_bin;
      distance_to_mid_res[i]=distance_to_mid_res[i]*pixels_per_bin;
      distance_to_mid[index]=-1;
    }
  free(x_field);
  free(y_field);
  free(radius_field);
  free(distance_to_mid);
  free(area_field);
  return;
}

double hux_heading(double delta_x, double delta_y)
{
  /*********************mHux_Heading*************************************
 - calculate: heading from change in x-position & y-position
 - requires: xchange (x), ychange (y)
 - returns: angle (0-359, 0=east, 90=north, 180=west, 270=south)
************************************************************************/
  double angle;
  if(delta_x==0)
    {
      if(delta_y==0) return(0);
      else if(delta_y>0) return(90);
      else return(270);
    }
  else if(delta_y==0)
    {
      if(delta_x<0) return(180);
      else return(0);
    }
  else 
    {
      angle = atanf(delta_y/delta_x)*57.29577951308232087685; /* 57.295... = 180/pi */
      if(delta_x<0) angle+=180;
      else if (delta_x>0&&delta_y<0) angle+=360;
      return(angle);
    }
}

void map_rotate(double* map,
		int x_bins,
		int y_bins,
		double deg,
		double invalid)
{
  //function to rotate a map around the middle of it
  int x;
  int y;
  double hypo;
  int x_adjacent=0;
  int y_opposite=0;
  double old_angle;
  double new_angle; // deg
  double new_angle_radian; //radians
  int new_index;
  int mid_x=x_bins/2;
  int mid_y=y_bins/2;
  int map_size=x_bins*y_bins;
  double* new_map;
  
  if (deg < 0 || deg > 360)
  {
    Rprintf("deg needs to be between 0 and 360 in map_rotate\n");
    return ;
  }
  new_map = (double*)malloc(map_size*sizeof(double));
  // set the new map to invalid
  set_array_to_value_double(new_map,map_size,invalid);
  
  // do the rotation
  for (int i = 0; i < map_size; i++)
  {
    get_x_and_y_bin_from_index(x_bins,
                               y_bins,
                               i,
                               &x,
                               &y);
                               
                               
    // relative to center
    x=x-mid_x;
    y=y-mid_y;
    hypo=sqrt(x*x+y*y);
    old_angle=hux_heading(x,y);
    new_angle=old_angle+deg;
    if (new_angle >= 360)
    {
      new_angle = new_angle -360;
    }
    new_angle_radian=degree_to_radian(new_angle);
    
    // now we need the angle from 0 to M_PI/2 (90 deg)
    while (new_angle_radian > M_PI/2)
    {
      new_angle_radian= new_angle_radian-M_PI/2;
    }
    
    // set the new coordinates for the pixel
    // add .5 and transform to int to round up to closest integer
    if (new_angle == 0)
    {
      x_adjacent=(int)(hypo+0.5);
      y_opposite=0;
    }
    if (new_angle > 0 && new_angle < 90)
    {
      x_adjacent=(int)((cos(new_angle_radian)*hypo)+.5);
      y_opposite=(int)((sin(new_angle_radian)*hypo)+.5);
    }
    if (new_angle == 90)
    {
      x_adjacent=0;
      y_opposite=(int)(hypo+.5);
    }
    if (new_angle > 90 && new_angle < 180)
    {
      x_adjacent=0-(int)((cos(M_PI/2-new_angle_radian)*hypo)+.5);
      y_opposite=(int)((sin(M_PI/2-new_angle_radian)*hypo)+.5);
    }
    if (new_angle == 180)
    {
      x_adjacent=0-(int)(hypo+.5);
      y_opposite=0;
    }
    if (new_angle > 180 && new_angle < 270)
    {
      x_adjacent=0-(int)((cos(new_angle_radian)*hypo)+.5);
      y_opposite=0-(int)((sin(new_angle_radian)*hypo)+.5);
    }
    if (new_angle == 270)
    {
      x_adjacent=0;
      y_opposite=0-(int)(hypo+.5);
    }
    if (new_angle > 270 && new_angle < 360)
    {
      x_adjacent=(int)((cos(M_PI/2-new_angle_radian)*hypo)+.5);
      y_opposite=0-(int)((sin(M_PI/2-new_angle_radian)*hypo)+.5);
    }
    
    // add the value of the center to the rotated coordinates
    x_adjacent=x_adjacent+mid_x;
    y_opposite=y_opposite+mid_y;
    
    // get the index of the bin in the array
    if (x_adjacent >= 0 && x_adjacent < x_bins && y_opposite>=0 && y_opposite < y_bins)
    {
      new_index=get_index_from_x_and_y_bin(x_bins,
                                           y_bins,
                                           x_adjacent,
                                           y_opposite);
      new_map[new_index]=map[i];
    }
    
  }
  for (int i = 0; i < map_size; i++)
  { map[i]=new_map[i];}
  free(new_map);
  return;
}


SEXP grid_score_cwrap(SEXP cell_lines_r,
			  SEXP auto_maps_r,
			  SEXP auto_num_bins_x_r,
			  SEXP auto_num_bins_y_r,
			  SEXP pixels_per_bin_r, 
			  SEXP number_fields_to_detect_r, // min number of bins in fields
			  SEXP min_num_bins_per_field_r,
			  SEXP field_threshold_r, 
			  SEXP invalid_r)
{
  int cell_lines = INTEGER_VALUE(cell_lines_r);
  double* all_autos = REAL(auto_maps_r);
  double* one_auto;
  int auto_num_bins_x = INTEGER_VALUE(auto_num_bins_x_r);
  int auto_num_bins_y = INTEGER_VALUE(auto_num_bins_y_r);
  int total_bins_auto= auto_num_bins_x*auto_num_bins_y;


  // create a copy of autocorrelation otherwise it will remove the fields from autocorrelation
  double* all_autos_copy = (double*)malloc(total_bins_auto*cell_lines*sizeof(double));
  double* one_auto_copy;
  
  
  for(int i = 0; i < cell_lines; i++){
    one_auto = all_autos + (i*total_bins_auto);
    one_auto_copy = all_autos_copy + (i*total_bins_auto);
    for(int x =0; x < auto_num_bins_x;x++)
      for(int y =0; y < auto_num_bins_y;y++)
        one_auto_copy[x*auto_num_bins_y+y]=one_auto[x+auto_num_bins_x*y]; // needs to be transpose for c code
  }
  
  SEXP out = PROTECT(allocVector(REALSXP,cell_lines));
  double* o = REAL(out);

  
  for(int i = 0; i < cell_lines; i++){
    one_auto=all_autos_copy+(i*total_bins_auto);
    o[i]= gridness_score(one_auto,
			    auto_num_bins_x,
			    auto_num_bins_y,
			    REAL(pixels_per_bin_r)[0],
			    INTEGER_VALUE(number_fields_to_detect_r),
			    INTEGER_VALUE(min_num_bins_per_field_r),
			    REAL(field_threshold_r)[0], 
			    REAL(invalid_r)[0]);
  }
  free(all_autos_copy);
  UNPROTECT(1);
  return (out);
}
double gridness_score(double* one_auto_map,
		      int auto_num_bins_x,
		      int auto_num_bins_y,
		      double pixels_per_bin,
		      int number_fields_to_detect, // min number of bins in fields
		      int min_num_bins_per_field,
		      double field_threshold, 
		      double invalid)
{
  // function to calculate the gridness score
  int num_fields = 7;
  int num_valid_fields;
  int num_bins = auto_num_bins_x*auto_num_bins_y;
  double* original_map;
  double* rotated_map;
  double* x_field;
  double* y_field;
  double* radius;
  double* distance_to_mid;
  int* area_field;
  int mid_x=auto_num_bins_x/2;
  int mid_y=auto_num_bins_y/2;
  int x;
  int y;
  int d;
  double min_radius; // for the ring
  double max_radius; // for the ring
  double corr_30;
  double corr_60;
  double corr_90;
  double corr_120;
  double corr_150;
  original_map = (double*)malloc(num_bins*sizeof(double));
  rotated_map =  (double*)malloc(num_bins*sizeof(double));
  x_field =  (double*)malloc(num_fields*sizeof(double));
  y_field =  (double*)malloc(num_fields*sizeof(double));
  radius =  (double*)malloc(num_fields*sizeof(double));
  distance_to_mid =  (double*)malloc(num_fields*sizeof(double));
  area_field =  (int*)malloc(num_fields*sizeof(int));
  // copy the original map before detection of fields
  for (int i = 0 ; i < num_bins ; i++)
    {
      original_map[i] = one_auto_map [i];
    }
  // get the position of the 6 fields close to mid of map
  grid_closest_peaks_to_middle(one_auto_map,
			       auto_num_bins_x,
			       auto_num_bins_y,
			       pixels_per_bin,
			       number_fields_to_detect,
			       min_num_bins_per_field,
			       field_threshold,
			       invalid,
			       x_field, // res for results
			       y_field,
			       radius,
			       distance_to_mid,
			       area_field,
			       num_fields);
  
  num_valid_fields=0;
  for(int i =0; i < num_fields; i++)
    {
      if(x_field[i]!=-1)
	{
	  num_valid_fields++;
	}
    }

  // to make the ring that will be rotated
  min_radius=radius[0]/pixels_per_bin; 
  if (min_radius >= auto_num_bins_x/4 || min_radius >= auto_num_bins_y/4) // this check needs to be relative to the largest distance of valid bin from center... see below.
    {
      min_radius=min_radius/2;
    }
  max_radius=0;
  for (int i = 1 ; i < num_fields; i++)
    {
      if (x_field[i] != -1)
	{
	  if ((distance_to_mid[i]+radius[i])/pixels_per_bin > max_radius)
	    {
	      max_radius=(distance_to_mid[i]+radius[i])/pixels_per_bin;
	    }
	}
    }

  if(max_radius==0)
    {
      max_radius=auto_num_bins_x/4;
    }
  
  if(max_radius<=min_radius)
    {
     // Rprintf("max_radius: %lf, min_radius: %lf, min_radius adjusted to %lf\n", max_radius, min_radius, max_radius/2);
      min_radius=max_radius/2;
    }

  // set the value outside the ring to invalid
  for (int i = 0 ; i < num_bins ; i++)
    {
      if (original_map[i] != invalid)
	{
	  // check the distance between that bin and the center of the map
	  get_x_and_y_bin_from_index(auto_num_bins_x,
				     auto_num_bins_y,
				     i,
				     &x,
				     &y);
	  d=(int)distance(mid_x,mid_y,x,y);
	  if (d<min_radius|| d>max_radius)
	    {
	      original_map[i] = invalid;
	    }
	}
    }
  
  // do the 5 correlations for the gridness score
  for (int j = 0; j < num_bins; j++)
    {
      rotated_map[j]=original_map[j];
    }
  map_rotate(rotated_map,
	     auto_num_bins_x,
	     auto_num_bins_y,
	     30,
	     invalid);
  corr_30=correlation(original_map,rotated_map,num_bins,invalid);
  for (int j = 0; j < num_bins; j++)
    {
      rotated_map[j]=original_map[j];
    }
  map_rotate(rotated_map,
	     auto_num_bins_x,
	     auto_num_bins_y,
	     60,
	     invalid);
  corr_60=correlation (original_map,rotated_map,num_bins,invalid);
  for (int j = 0; j < num_bins; j++)
    {
      rotated_map[j]=original_map[j];
    }
  map_rotate(rotated_map,
	     auto_num_bins_x,
	     auto_num_bins_y,
	     90,
	     invalid);
  corr_90=correlation (original_map,rotated_map,num_bins,invalid);
  for (int j = 0; j < num_bins; j++)
    {
      rotated_map[j]=original_map[j];
    }
  map_rotate(rotated_map,
	     auto_num_bins_x,
	     auto_num_bins_y,
	     120,
	     invalid);
  corr_120=correlation (original_map,rotated_map,num_bins,invalid);
  for (int j = 0; j < num_bins; j++)
    {
      rotated_map[j]=original_map[j];
    }
  map_rotate(rotated_map,
	     auto_num_bins_x,
	     auto_num_bins_y,
	     150,
	     invalid);
  corr_150=correlation (original_map,rotated_map,num_bins,invalid);
  
  free(original_map);
  free(rotated_map);
  free(x_field);
  free(y_field);
  free(radius);
  free(distance_to_mid);
  free(area_field);
  return ((corr_60+corr_120)/2) - ((corr_30+corr_90+corr_150)/3);
}

SEXP grid_orientation_cwrap(SEXP cell_lines_r,
			    SEXP auto_maps_r,
			    SEXP auto_num_bins_x_r,
			    SEXP auto_num_bins_y_r,
			    SEXP pixels_per_bin_r,
			    SEXP number_fields_to_detect_r,
			    SEXP min_num_bins_per_field_r,
			    SEXP field_threshold_r,
			    SEXP invalid_r)
{
  int cell_lines = INTEGER_VALUE(cell_lines_r);
  double* all_autos = REAL(auto_maps_r);
  double* one_auto;
  int auto_num_bins_x = INTEGER_VALUE(auto_num_bins_x_r);
  int auto_num_bins_y = INTEGER_VALUE(auto_num_bins_y_r);
  int total_bins_auto= auto_num_bins_x*auto_num_bins_y;


  // create a copy of autocorrelation otherwise it will remove the fields from autocorrelation
  double* all_autos_copy = (double*)malloc(total_bins_auto*cell_lines*sizeof(double));
  for(int i =0; i <total_bins_auto*cell_lines;i++)
    all_autos_copy[i]=all_autos[i];

  SEXP out = PROTECT(allocVector(REALSXP,cell_lines));
  double* o = REAL(out);
  
  for(int i = 0; i < cell_lines; i++){
    one_auto=all_autos_copy+(i*total_bins_auto);
    o[i]= grid_orientation(one_auto,
			   auto_num_bins_x,
			   auto_num_bins_y,
			   REAL(pixels_per_bin_r)[0],
			   INTEGER_VALUE(number_fields_to_detect_r),
			   INTEGER_VALUE(min_num_bins_per_field_r),
			   REAL(field_threshold_r)[0], 
			   REAL(invalid_r)[0]);
  }
  free(all_autos_copy);
  UNPROTECT(1);
  return (out);

}



double grid_orientation(double *one_auto,
			int x_bins,
			int y_bins,
			double pixels_per_bin,
			int num_fields_to_detect,
			int min_num_bins_fields,
			float threshold,
			double invalid)
{
  double* x_field;
  double* y_field;
  double* radius;
  double* distance_to_mid;
  int* area_field;
  double* orientation;
  int num_fields = 7;
  double smallest_orientation=0;
  int mid_x=(int)(x_bins*pixels_per_bin/2);
  int mid_y=(int)(y_bins*pixels_per_bin/2);

  x_field = (double*)malloc(num_fields*sizeof(double));
  y_field = (double*)malloc(num_fields*sizeof(double));
  radius = (double*)malloc(num_fields*sizeof(double));
  distance_to_mid =(double*)malloc(num_fields*sizeof(double));
  area_field = (int*)malloc(num_fields*sizeof(int));
  orientation = (double*)malloc(num_fields*sizeof(double));

  grid_closest_peaks_to_middle(one_auto,
			       x_bins,
			       y_bins,
			       pixels_per_bin,
			       num_fields_to_detect,
			       min_num_bins_fields,
			       threshold,
			       invalid,
			       x_field, // res for results
			       y_field,
			       radius,
			       distance_to_mid,
			       area_field,
			       num_fields);

  // we now calculate the angle between center of the map and each field
  smallest_orientation=360;
  for (int i = 1; i < num_fields; i++) // start at one to eliminate the center of map
    {
      if (x_field[i]!=-1)
	{
	  // now we need to get the angle for the 
	  orientation[i]=hux_heading(x_field[i]-mid_x,y_field[i]-mid_y);
	  if (orientation[i]<smallest_orientation)
	    {
	      smallest_orientation=orientation[i];
	    }
	}
    }
  free(x_field);
  free(y_field);
  free(radius);
  free(distance_to_mid);
  free(area_field);
  free(orientation);
  return smallest_orientation;
}


SEXP grid_spacing_cwrap(SEXP cell_lines_r,
			    SEXP auto_maps_r,
			    SEXP auto_num_bins_x_r,
			    SEXP auto_num_bins_y_r,
			    SEXP pixels_per_bin_r,
			    SEXP number_fields_to_detect_r,
			    SEXP min_num_bins_per_field_r,
			    SEXP field_threshold_r,
			    SEXP invalid_r)
{
  int cell_lines = INTEGER_VALUE(cell_lines_r);
  double* all_autos = REAL(auto_maps_r);
  double* one_auto;
  int auto_num_bins_x = INTEGER_VALUE(auto_num_bins_x_r);
  int auto_num_bins_y = INTEGER_VALUE(auto_num_bins_y_r);
  int total_bins_auto= auto_num_bins_x*auto_num_bins_y;
  
  // create a copy of autocorrelation to transpose and to avoid removing the fields from original autocorrelation
  double* all_autos_copy = (double*)malloc(total_bins_auto*cell_lines*sizeof(double));
  double* one_auto_copy;
  
  for(int i = 0; i < cell_lines; i++){
    one_auto = all_autos + (i*total_bins_auto);
    one_auto_copy = all_autos_copy + (i*total_bins_auto);
    for(int x =0; x < auto_num_bins_x;x++)
      for(int y =0; y < auto_num_bins_y;y++)
        one_auto_copy[x*auto_num_bins_y+y]=one_auto[x+auto_num_bins_x*y]; // needs to be transpose for c code
  }
  // vector of spacing to return to R
  SEXP out = PROTECT(allocVector(REALSXP,cell_lines));
  double* o = REAL(out);
  
  for(int i = 0; i < cell_lines; i++){
    one_auto=all_autos_copy+(i*total_bins_auto);
    o[i]= grid_spacing(one_auto,
		       auto_num_bins_x,
		       auto_num_bins_y,
		       REAL(pixels_per_bin_r)[0],
		       INTEGER_VALUE(number_fields_to_detect_r),
		       INTEGER_VALUE(min_num_bins_per_field_r),
		       REAL(field_threshold_r)[0], 
		       REAL(invalid_r)[0]);
  }
  
 
  free(all_autos_copy);
  UNPROTECT(1);
  return (out);
}

double grid_spacing(double *map,
		    int x_bins,
		    int y_bins,
		    double pixels_per_bin,
		    int num_fields_to_detect,
		    int min_num_bins_fields,
		    float threshold,
		    double invalid)
{
  /********************************************************
      funtion to detect fields in a map and remove them
  **********************************************************/
  
  double* x_field;
  double* y_field;
  double* radius;
  double* distance_to_mid;
  int* area_field;
  double mean_distance_to_mid =0;
  int num_fields = 7;
  int valid_fields=0;

  x_field = (double*)malloc(num_fields*sizeof(double));
  y_field = (double*)malloc(num_fields*sizeof(double));
  radius = (double*)malloc(num_fields*sizeof(double));
  distance_to_mid =(double*)malloc(num_fields*sizeof(double));
  area_field = (int*)malloc(num_fields*sizeof(int));

  grid_closest_peaks_to_middle(map,
			       x_bins,
			       y_bins,
			       pixels_per_bin,
			       num_fields_to_detect,
			       min_num_bins_fields,
			       threshold,
			       invalid,
			       x_field, // res for results
			       y_field,
			       radius,
			       distance_to_mid,
			       area_field,
			       num_fields);
  
  for (int i = 1; i < num_fields; i++) // start to 1 to eliminate the field at center
    {
      if (distance_to_mid[i] != -1)
	{
	  mean_distance_to_mid=mean_distance_to_mid+distance_to_mid[i];
	  valid_fields++;
	}
    }
  
  free(x_field);
  free(y_field);
  free(radius);
  free(distance_to_mid);
  free(area_field);
  return mean_distance_to_mid/valid_fields;
}

void remove_internal_bins_from_border(int num_bins_x, int num_bins_y, int* border_map, int* border_x, int* border_y, int* num_bins_border)
{
  // takes care of the "border" bins that are in the middle of the map
  
  // is_border needs to be above 0 for a border pixel to be kept
  int min;
  int max;
  int min_index;
  int max_index;
  int* is_border= (int*) malloc(sizeof(int)* (*num_bins_border));
  for(int i = 0; i < (*num_bins_border); i++)
    is_border[i]=0;

  // find the outside bins for each row
  for(int x = 0; x < num_bins_x; x++)
  {
    min=num_bins_x;
    max=0;
    min_index=-1;
    max_index=-1;
    for(int i = 0; i < (*num_bins_border);i++)
    {
      if(border_x[i]==x)
      {
        if(border_y[i]>max){
          max=border_y[i];
          max_index=i;
        }
        if(border_y[i]<min){
          min=border_y[i];
          min_index=i;
        }
      }
    }
    // for this row, keep only the extrem values
    if(max_index!=-1)
      is_border[max_index]=1;
    if(min_index!=-1)
      is_border[min_index]=1;
  }
  // find the outside bins for each column
  for(int y = 0; y < num_bins_y; y++)
  {
    min=num_bins_y;
    max=0;
    min_index=-1;
    max_index=-1;
    for(int i = 0; i < (*num_bins_border);i++)
    {
      if(border_y[i]==y)
      {
        if(border_x[i]>max){
          max=border_x[i];
          max_index=i;
        }
        if(border_x[i]<min){
          min=border_x[i];
          min_index=i;
        }
      }
    }
     // for this row, keep only the extrem values
     if(max_index!=-1)
       is_border[max_index]=1;
     if(min_index!=-1)
       is_border[min_index]=1;
    }

    int new_num_bins=0;
  for(int i = 0; i < (*num_bins_border);i++)
  {
    if(is_border[i]>0)
      new_num_bins++;
  }

  
  // remove the bins that is_border==0
  int removed=0;
  for(int i = 0; i < (*num_bins_border);i++)
  {
    if(is_border[i]==0)
    {
      // remove from the map
      border_map[(border_x[i]*num_bins_y)+border_y[i]]=0;
      removed++;
    }
    else
    {
      border_x[i-removed]=border_x[i];
      border_y[i-removed]=border_y[i];
    }
  }
  *num_bins_border=new_num_bins;
  free(is_border);
}


int identify_border_pixels_in_occupancy_map(double* occ_map, int num_bins_x, int num_bins_y,int* border_map, int* border_x, int* border_y, int* num_bins_border)
{
  for(int i = 0; i < num_bins_x*num_bins_y;i++) // set border map to 0
    {border_map[i]=0;}
  *num_bins_border=0;
  while(find_border_starting_point(occ_map, num_bins_x, num_bins_y,border_map,border_x,border_y,num_bins_border))
    {
      //recursive algorhythm! Inspired by conversation with Jozsef and Catherine in Vienna
      find_an_adjacent_border_pixel(occ_map, num_bins_x, num_bins_y,border_map,border_x,border_y,num_bins_border);
    }
  
  // remove internal bins 
  remove_internal_bins_from_border(num_bins_x, num_bins_y, border_map, border_x, border_y, num_bins_border);
  return 0;
}
int find_border_starting_point(double* occ_map, int num_bins_x, int num_bins_y,int*border_map,int*border_x,int* border_y,int* num_bins_border)
{
  for (int x = 0 ; x < num_bins_x; x++)
  {
    for (int y = 0; y < num_bins_y; y++)
    {
      if (x!=0&&x!=num_bins_x-1&&y!=0&&y!=num_bins_y-1)
      {  
        // a starting point could be a valid pixel, with 3 non visited pixels on right, and none on the left, or reverse
        if((border_map[((x)*num_bins_y)+y]!=1 &&
           occ_map[((x)*num_bins_y)+y]!=-1 &&
           occ_map[((x-1)*num_bins_y)+y-1]==-1 &&
           occ_map[((x-1)*num_bins_y)+y]==-1 &&
           occ_map[((x-1)*num_bins_y)+y+1]==-1 &&
           occ_map[((x+1)*num_bins_y)+y-1]!=-1 &&
           occ_map[((x+1)*num_bins_y)+y]!=-1 && 
           occ_map[((x+1)*num_bins_y)+y+1]!=-1)
             ||
               (border_map[((x)*num_bins_y)+y]!=1 &&
               occ_map[((x)*num_bins_y)+y]!=-1 &&
               occ_map[((x-1)*num_bins_y)+y-1]!=-1 &&
               occ_map[((x-1)*num_bins_y)+y]!=-1 &&
               occ_map[((x-1)*num_bins_y)+y+1]!=-1 &&
               occ_map[((x+1)*num_bins_y)+y-1]==-1 &&
               occ_map[((x+1)*num_bins_y)+y]==-1 &&
               occ_map[((x+1)*num_bins_y)+y+1]==-1))
        {
          border_map[(x*num_bins_y)+y]=1;
          border_x[*num_bins_border]=x;
          border_y[*num_bins_border]=y;
         // Rprintf("starting at %d %d\n",x,y);
          (*num_bins_border)++;
          return 1;
        }
        // a border pixels could have 3 non visited pixels above, and none below, or reverse
        if((border_map[((x)*num_bins_y)+y]!=1 &&
           occ_map[((x)*num_bins_y)+y]!=-1 &&
           occ_map[((x-1)*num_bins_y)+y-1]==-1 &&
           occ_map[((x)*num_bins_y)+y-1]==-1 &&
           occ_map[((x+1)*num_bins_y)+y-1]==-1 &&
           occ_map[((x-1)*num_bins_y)+y+1]!=-1 &&
           occ_map[((x)*num_bins_y)+y+1]!=-1 &&
           occ_map[((x+1)*num_bins_y)+y+1]!=-1)
             ||
               (border_map[((x)*num_bins_y)+y]!=1 &&
               occ_map[((x)*num_bins_y)+y]!=-1 &&
               occ_map[((x-1)*num_bins_y)+y-1]!=-1 &&
               occ_map[((x)*num_bins_y)+y-1]!=-1 &&
               occ_map[((x+1)*num_bins_y)+y-1]!=-1 &&
               occ_map[((x-1)*num_bins_y)+y+1]==-1 &&
               occ_map[((x)*num_bins_y)+y+1]==-1 &&
               occ_map[((x+1)*num_bins_y)+y+1]==-1))
        {
          border_map[(x*num_bins_y)+y]=1;
          border_x[*num_bins_border]=x;
          border_y[*num_bins_border]=y;
          (*num_bins_border)++;
         // Rprintf("starting at %d %d\n",x,y);
          return 1;
        }
        
        // a border pixel could have 3 non visited pixels and be at one of the 4 corners of the rectangle
        // need this otherwise we miss 2 corners in a rectangular map
        if((border_map[((x)*num_bins_y)+y]!=1 && // top left
           occ_map[((x)*num_bins_y)+y]!=-1 &&
           occ_map[((x-1)*num_bins_y)+y]==-1 &&
           occ_map[((x-1)*num_bins_y)+y+1]==-1 &&
           occ_map[((x)*num_bins_y)+y+1]==-1 &&
           occ_map[((x)*num_bins_y)+y-1]!=-1 &&
           occ_map[((x+1)*num_bins_y)+y]!=-1 &&
           occ_map[((x+1)*num_bins_y)+y+1]!=-1)
          ||
            (border_map[((x)*num_bins_y)+y]!=1 && // top right
            occ_map[((x)*num_bins_y)+y]!=-1 &&
            occ_map[((x)*num_bins_y)+y+1]==-1 &&
            occ_map[((x+1)*num_bins_y)+y+1]==-1 &&
            occ_map[((x+1)*num_bins_y)+y]==-1 &&
            occ_map[((x-1)*num_bins_y)+y]!=-1 &&
            occ_map[((x-1)*num_bins_y)+y-1]!=-1 &&
            occ_map[((x)*num_bins_y)+y-1]!=-1)
          ||
            (border_map[((x)*num_bins_y)+y]!=1 && // bottom left
            occ_map[((x)*num_bins_y)+y]!=-1 &&
            occ_map[((x-1)*num_bins_y)+y]==-1 &&
            occ_map[((x-1)*num_bins_y)+y-1]==-1 &&
            occ_map[((x)*num_bins_y)+y-1]==-1 &&
            occ_map[((x)*num_bins_y)+y+1]!=-1 &&
            occ_map[((x+1)*num_bins_y)+y+1]!=-1 &&
            occ_map[((x+1)*num_bins_y)+y]!=-1)
            ||
            (border_map[((x)*num_bins_y)+y]!=1 && // bottom right
            occ_map[((x)*num_bins_y)+y]!=-1 &&
            occ_map[((x)*num_bins_y)+y-1]==-1 &&
            occ_map[((x+1)*num_bins_y)+y-1]==-1 &&
            occ_map[((x+1)*num_bins_y)+y]==-1 &&
            occ_map[((x-1)*num_bins_y)+y]!=-1 &&
            occ_map[((x-1)*num_bins_y)+y+1]!=-1 &&
            occ_map[((x)*num_bins_y)+y+1]!=-1))
        {
          border_map[(x*num_bins_y)+y]=1;
          border_x[*num_bins_border]=x;
          border_y[*num_bins_border]=y;
          (*num_bins_border)++;
         // Rprintf("starting at %d %d\n",x,y);
          return 1;
        }
      }
      else
      {// border could be a visited bin on the edge of the occ map
        if (occ_map[(x*num_bins_y)+y]!=-1&&border_map[((x)*num_bins_y)+y]!=1)
        {
          border_map[(x*num_bins_y)+y]=1;
          border_x[*num_bins_border]=x;
          border_y[*num_bins_border]=y;
          (*num_bins_border)++;
        //  Rprintf("starting at %d %d\n",x,y);
          return 1;
        }
      }
    }
  }
  return 0;
}

int find_an_adjacent_border_pixel(double* occ_map, int num_bins_x, int num_bins_y,int*border_map,int*border_x,int* border_y,int* num_bins_border)
{
  // look for a pixels around the last added pixel that could be a border, in the 9 pixels around it
  for (int x = border_x[(*num_bins_border)-1]-1;x<=border_x[(*num_bins_border)-1]+1;x++)
  {
    for (int y = border_y[(*num_bins_border)-1]-1;y<=border_y[(*num_bins_border)-1]+1;y++)
    {
      if(
        ((x>0 && // not in x==0, we need to see a -1.0 value on the left side
          x<num_bins_x-1 && // not is x = length(x), we need to see a -1.0 value on the right side
          border_map[((x)*num_bins_y)+y]!=1 &&  // not already considered a border pixel
          occ_map[((x)*num_bins_y)+y]!=-1&&y<num_bins_y) && // bin was visited by animal and is not after y limit
          (occ_map[((x-1)*num_bins_y)+y]==-1 || occ_map[((x+1)*num_bins_y)+y]==-1)) // should have an invalid bin on right or left side
        || 
          ((y>0 && 
          y<num_bins_y-1 && 
          border_map[((x)*num_bins_y)+y]!=1 && 
          occ_map[((x)*num_bins_y)+y]!=-1&&x<num_bins_x) &&
          (occ_map[(x*num_bins_y)+y-1]==-1 || occ_map[(x*num_bins_y)+y+1]==-1)))
      {
        border_map[(x*num_bins_y)+y]=1;
        border_x[*num_bins_border]=x;
        border_y[*num_bins_border]=y;
     //   Rprintf("adding %d %d\n",x,y);
        (*num_bins_border)++;
        find_an_adjacent_border_pixel(occ_map,num_bins_x,num_bins_y,border_map,border_x,border_y,num_bins_border);// recursive search
        return 1;
      }
    }
  }
  return 0;
}


int assign_wall_to_border_pixels(int num_bins_x, int num_bins_y, int* border_x, int* border_y, int* num_bins_border,int* wall_id,int* border_map)
{
  /* 
     function assumes that there are 2 vertical and 2 horizontal walls.
     the two vertical or horizontal walls should be in different half of the map
     We simply count the number of border pixels for each row or column of the map
     Walls should results in a high number of pixel for a given row or column
   
     Note that a pixel can only be assigned to one border. It is assigned to the closest wall, if distances are equal the assignation
     is arbitrary set by the order of if statements
  */
  
  double* col_sum = (double*) malloc(num_bins_x*sizeof(double));
  double* row_sum = (double*) malloc(num_bins_y*sizeof(double));
  double sum,max;
  int h1=0,h2=0,v1=0,v2=0; // coordinate of the horizontal and vertical walls
  int dist_h1=0;
  int dist_h2=0;
  int dist_v1=0;
  int dist_v2=0;
  int* x;
  int* y;
  int* id;
  int num_bins_wall;

  for(int i = 0; i < num_bins_x; i++)
  {
    sum=0;
    for(int j = 0; j < *num_bins_border; j++)
    {
      if(border_x[j]==i)
        sum++;
    }
    col_sum[i]=sum;
  }
  
  for(int i = 0; i < num_bins_y; i++)
  {
    sum=0;
    for(int j = 0; j < *num_bins_border; j++)
    {
      if(border_y[j]==i)
        sum++;
    }
    row_sum[i]=sum;
  }
  
  max=0;
  for(int i = 0; i < num_bins_x/2; i++)
  {
    if(col_sum[i]>max)
    {
      max=col_sum[i];
      v1=i;
    }
  }
  max=0;
  for(int i = num_bins_x/2; i < num_bins_x; i++)
  {
    if(col_sum[i]>max)
    {
      max=col_sum[i];
      v2=i;
    }
  }
  max=0;
  for(int i = 0; i < num_bins_y/2; i++)
  {
    if(row_sum[i]>max)
    {
      max=row_sum[i];
      h1=i;
    }
  }
  max=0;
  for(int i = num_bins_y/2; i < num_bins_y; i++)
  {
    if(row_sum[i]>max)
    {
      max=row_sum[i];
      h2=i;
    }
  }
  
  for(int i =0; i < *num_bins_border;i++)
    wall_id[i]=-1;
  
  for(int i =0; i < *num_bins_border;i++)
  {
    dist_v1=sqrt((border_x[i]-v1)*(border_x[i]-v1));
    dist_v2=sqrt((border_x[i]-v2)*(border_x[i]-v2));
    dist_h1=sqrt((border_y[i]-h1)*(border_y[i]-h1));
    dist_h2=sqrt((border_y[i]-h2)*(border_y[i]-h2));
 
 
    // assign to the closest border
    if(dist_v1<2&&dist_v1<=dist_v2&&dist_v1<=dist_h1&&dist_v1<=dist_h2)
       wall_id[i]=0;
    if(dist_v2<2&&dist_v2<=dist_v1&&dist_v2<=dist_h1&&dist_v2<=dist_h2)
      wall_id[i]=1;
    if(dist_h1<2&&dist_h1<=dist_h2&&dist_h1<=dist_v1&&dist_h1<=dist_v2)
      wall_id[i]=2;
    if(dist_h2<2&&dist_h2<=dist_h1&&dist_h2<=dist_v1&&dist_h2<=dist_v2)
      wall_id[i]=3;
  }
  
  // if a pixels is not associated to a wall, remove it from the border_map;
  for(int i =0; i < *num_bins_border;i++)
  {
    if(wall_id[i]==-1)
    {
      border_map[(border_x[i]*num_bins_y)+border_y[i]]=0;
    }
  }
  // remove border pixels that are not close to a wall
  x = (int*)malloc((*num_bins_border)*sizeof(int));
  y = (int*)malloc((*num_bins_border)*sizeof(int));
  id =(int*)malloc((*num_bins_border)*sizeof(int));
  num_bins_wall=0;
  
  for(int i =0; i < *num_bins_border;i++)
  {
    if(wall_id[i]!=-1)
    {
      x[num_bins_wall]=border_x[i];
      y[num_bins_wall]=border_y[i];
      id[num_bins_wall]=wall_id[i];
      num_bins_wall++;
    }
  }
  for(int i = 0; i < num_bins_wall; i++)
  {
    border_x[i]=x[i];
    border_y[i]=y[i];
    wall_id[i]=id[i];
  }
  
  (*num_bins_border)=num_bins_wall;
  
  free(col_sum);
  free(row_sum);
  free(x);
  free(y);
  free(id);
  return 0;
}




SEXP border_detection_rectangular_environment_cwrap(SEXP num_bins_x_r,
                                                SEXP num_bins_y_r,
                                                SEXP occ_map_r)
{
  int num_bins_x = INTEGER_VALUE(num_bins_x_r);
  int num_bins_y = INTEGER_VALUE(num_bins_y_r);
  double* occ_map = REAL(occ_map_r);
  
  int total_bins= num_bins_x*num_bins_y;
  
  // transpose the occ map for c functions
  double* tocc_map = (double*)malloc(total_bins*sizeof(double));
  for(int x = 0; x < num_bins_x; x++) 
    for(int y = 0; y < num_bins_y; y++)
      tocc_map[x*num_bins_y+y]=occ_map[x+num_bins_x*y];
  
  // 1) identify the borders in the occupancy map
  int num_bins_border=0;
  SEXP out = PROTECT(allocMatrix(INTSXP,num_bins_x,num_bins_y));
  int* border_map = (int*) malloc(total_bins*sizeof(int));
  int* border_x =  (int*) malloc(total_bins*sizeof(int)); // list of x index for border bins
  int* border_y = (int*) malloc(total_bins*sizeof(int)); // list of y index for border bins
 // Rprintf("identify_border_pixels_in_occupancy_map() call\n");
  identify_border_pixels_in_occupancy_map(tocc_map,num_bins_x,num_bins_y,border_map,border_x,border_y,&num_bins_border);
  
  // 2) identify the four walls, remove border pixels that are not next to one of the 4 walls
  int* wall_id =  (int*) malloc(num_bins_border*sizeof(int));
  //Rprintf("assign_wall_to_border_pixel() call\n");
  assign_wall_to_border_pixels(num_bins_x, num_bins_y, border_x, border_y, &num_bins_border, wall_id,border_map);
  
  int sum;
  for(int j=0; j<4;j++){
    sum=0;
    for(int i = 0; i < num_bins_border;i++)
    {
      if(wall_id[i]==j)
        sum++;      
    }
   // Rprintf("wall id: %d, sum: %d\n",j,sum);
  }
  
  for(int i = 0; i < num_bins_border;i++)
  {
    border_map[border_x[i]*num_bins_y+border_y[i]]=wall_id[i]+1;
  }
  
  int* o = INTEGER_POINTER(out);
  //transpose the maps for R
  for(int x =0; x < num_bins_x;x++)
    for(int y =0; y < num_bins_y;y++)
        o[y*num_bins_x+x]=border_map[x*num_bins_y+y];
  

  free(border_map);
  free(border_x);
  free(border_y);
  free(tocc_map);
  UNPROTECT(1);
  return (out);
  
}


SEXP border_score_rectangular_environment_cwrap(SEXP cells_r,
						SEXP cell_lines_r,
						SEXP num_bins_x_r,
						SEXP num_bins_y_r,
						SEXP occ_map_r,
						SEXP maps_r,
						SEXP percent_threshold_field_r,
						SEXP min_bins_in_field_r)
{
  int* cells = INTEGER_POINTER(cells_r);
  int cell_lines = INTEGER_VALUE(cell_lines_r);
  int num_bins_x = INTEGER_VALUE(num_bins_x_r);
  int num_bins_y = INTEGER_VALUE(num_bins_y_r);
  double* occ_map = REAL(occ_map_r);
  double* maps = REAL(maps_r);
  double* one_map;
  double* one_map_copy;

  int total_bins= num_bins_x*num_bins_y;

  // create a copy of map otherwise it will remove the fields from it
  double* maps_copy = (double*)malloc(total_bins*cell_lines*sizeof(double));
  double* tocc_map = (double*)malloc(total_bins*sizeof(double));
  // transpose the map
  for(int i =0; i < cell_lines;i++){
    one_map=maps+(i*total_bins);
    one_map_copy=maps_copy+(i*total_bins);
    for(int x = 0; x < num_bins_x; x++) 
      for(int y = 0; y < num_bins_y; y++)
      one_map_copy[x*num_bins_y+y]=one_map[x+num_bins_x*y]; 
  }
  for(int x = 0; x < num_bins_x; x++) 
    for(int y = 0; y < num_bins_y; y++)
      tocc_map[x*num_bins_y+y]=occ_map[x+num_bins_x*y];
  

  SEXP out = PROTECT(allocMatrix(REALSXP,4,cell_lines));
  double* o = REAL(out);

  double* cm = (double*) malloc(cell_lines*sizeof(double));
  double* dm = (double*) malloc(cell_lines*sizeof(double));
  double* border_score =  (double*) malloc(cell_lines*sizeof(double));
  int* num_fields_detected = (int*) malloc(cell_lines*sizeof(int));

  border_score_rectangular_environment(cells,
				       cell_lines,
				       num_bins_x,
				       num_bins_y,
				       tocc_map,
				       maps_copy,
				       REAL(percent_threshold_field_r)[0],
				       INTEGER_VALUE(min_bins_in_field_r),
				       border_score, // border score
				       cm,
				       dm,
				       num_fields_detected);
  
  for(int i = 0; i < cell_lines; i++)
    {
      o[i*4+0]=border_score[i];
      o[i*4+1]=cm[i];
      o[i*4+2]=dm[i];
      o[i*4+3]=(double)num_fields_detected[i];

    }

  free(maps_copy);
  free(tocc_map);
  free(cm);
  free(dm);
  free(num_fields_detected);

  UNPROTECT(1);
  return (out);

}



void border_score_rectangular_environment(int* cells,
					  int cell_lines,
					  int num_bins_x,
					  int num_bins_y,
					  double* occ_map,
					  double* maps,
					  double percent_threshold_field,
					  int min_bins_in_field,
					  double* border_score,
					  double* cm,
					  double* dm,
					  int* num_fields_detected)


{
  // 1) identify the borders in the occupancy map
  int num_bins_border=0;
  int total_bins=num_bins_x*num_bins_y;
  int* border_map = (int*) malloc(total_bins*sizeof(int)); // map
  int* border_x =  (int*) malloc(total_bins*sizeof(int)); // list of x index for border bins
  int* border_y = (int*) malloc(total_bins*sizeof(int)); // list of y index for border bins
  identify_border_pixels_in_occupancy_map(occ_map,num_bins_x,num_bins_y,border_map,border_x,border_y,&num_bins_border);
  
  // 2) identify the four walls, remove border pixels that are not next to one of the 4 walls
  int* wall_id =  (int*) malloc(num_bins_border*sizeof(int));
  assign_wall_to_border_pixels(num_bins_x, num_bins_y, border_x, border_y, &num_bins_border, wall_id,border_map);
  
  double* one_map;
  // for each valid pixel in occ_map, find the closest distance to a border pixel, get the largest closest distance
  double min_distance_one_pixel=0;
  double dist=0;
  int number_fields_detected=0;
  double* bin_distance_to_nearest_border = (double*)malloc(total_bins*sizeof(double));
  double CM=0;
  double max_CM=0;
  double DM=0;
  int num_pixels_one_wall;
  int pixels_covered_one_wall;
  double* one_field_map; // rate of one firing field
  double* all_fields_map; // rate of all firing fields
  double* detection_map; // where fields get set to -1
  double max_fr_remaining=0;
  double max_fr;
  double threshold_hz;
  double mean_x_field;
  double mean_y_field;
  double radius_field;
  int num_bins_field;
  double sum_firing_rate_in_fields;
  
  one_field_map =  (double*)malloc(total_bins*sizeof(double));
  all_fields_map = (double*)malloc(total_bins*sizeof(double));
  detection_map = (double*)malloc(total_bins*sizeof(double));
  
  double max_possible_distance_to_border=0;
  for(int x = 0; x < num_bins_x; x++)
  {
    for(int y = 0; y < num_bins_y; y++)
    {
      if (occ_map[(x*num_bins_y) + y] != -1)
      {
        min_distance_one_pixel=num_bins_x+num_bins_y;
        for (int i = 0; i < num_bins_border; i++)
        {
          dist=distance(x,y,border_x[i],border_y[i]);
          if(dist<min_distance_one_pixel)
          {
            min_distance_one_pixel=dist;
          }
        }
        if(min_distance_one_pixel>max_possible_distance_to_border)
        {
          max_possible_distance_to_border=min_distance_one_pixel;
        }
      }
    }
  }
  
  for(int i = 0; i < cell_lines; i++)
  { // for each cell, 
    
    max_CM=0;
    one_map = maps + (i*total_bins);	      
    for(int j = 0; j < total_bins; j++)
    {
      one_field_map[j] =-1.0;
      all_fields_map[j] =-1.0;
      detection_map[j]=one_map[j];
    }
    
    // calculate the maximum CM
    max_fr=find_max_double(total_bins,detection_map);
    max_fr_remaining=max_fr;
    threshold_hz = max_fr*percent_threshold_field/100;
    number_fields_detected=0;
    //	      cerr << "max fr: " << max_fr << '\n';
    while(max_fr_remaining>=threshold_hz)
    {
      //		  cerr << "detect field with peak rate " << max_fr << '\n';
      detect_one_field_with_field(detection_map,
                                  num_bins_x,
                                  num_bins_y,
                                  min_bins_in_field,
                                  threshold_hz,
                                  &mean_x_field,
                                  &mean_y_field,
                                  &radius_field,
                                  &num_bins_field,
                                  -1.0,
                                  one_field_map);
      if (mean_x_field != -1)
      {// valid field detected
        
        number_fields_detected++;
        // copy the field to the all field maps
        for(int k=0; k < total_bins;k++)
        {
          if(one_field_map[k]!=-1.0)
            all_fields_map[k]=one_field_map[k];
        }
        
        // get the fraction of pixels along a wall that was occupied by the field.
        for(int k =0; k < 4; k++) // for the 4 walls
        {
          num_pixels_one_wall=0;
          pixels_covered_one_wall=0;
          // number of pixels for this wall
          for(int m =0; m < num_bins_border; m++)
          {
            if(wall_id[m]==k)
              num_pixels_one_wall++;
          }
          // number of occupied wall pixels
          for(int m =0; m < num_bins_border; m++)
          {
            if(wall_id[m]==k)
              if(one_field_map[(border_x[m]*num_bins_y)+border_y[m]]!=-1.0)
                pixels_covered_one_wall++;
          }
          
          if(num_pixels_one_wall!=0)
            CM=(double)pixels_covered_one_wall/(double)num_pixels_one_wall;
          else
            CM=0;
          
          if(CM>max_CM)
            max_CM=CM;
          //			  cerr << "wall " << k << " num_pix:" << num_pixels_one_wall << " cov_pix:" << pixels_covered_one_wall << " CM:" << CM << " fieldsize:" << num_bins_field << '\n'; 
        }
      }
      max_fr_remaining=find_max_double(total_bins,detection_map);
    }
    CM=max_CM;
    
    // calculate DM
    sum_firing_rate_in_fields=sum_double(total_bins,all_fields_map,-1.0);
    DM=0;
    for(int x = 0; x < num_bins_x; x++)
    {
      for(int y = 0; y < num_bins_y; y++)
      {
        if (all_fields_map[(x*num_bins_y) + y] != -1.0)// field bin
        {
          // find the closest distance to a border
          min_distance_one_pixel=num_bins_x+num_bins_y;
          for (int j = 0; j < num_bins_border; j++)
          {
            dist=distance(x,y,border_x[j],border_y[j]);
            if(dist<min_distance_one_pixel)
            {
              min_distance_one_pixel=dist;
            }
          }
          //  cerr << min_distance_one_pixel << " " << max_possible_distance_to_border << " " << one_place_map[(x*num_bins_y)+y] << '\n';
          DM+=(min_distance_one_pixel/max_possible_distance_to_border)*one_map[(x*num_bins_y) + y];
        }
      }
    }
    DM=DM/sum_firing_rate_in_fields;
    if(number_fields_detected>0)
    {
      cm[i]= CM;
      dm[i]= DM;
      num_fields_detected[i]=number_fields_detected;
      border_score[i]=(CM-DM)/(CM+DM);
    }
    else
    {
      cm[i]= NAN;
      dm[i]= NAN;
      num_fields_detected[i]=number_fields_detected;
      border_score[i]=(CM-DM)/(CM+DM);
    }
  }
  free(bin_distance_to_nearest_border);
  free(one_field_map);
  free(all_fields_map);
  free(detection_map);
  free(wall_id);
}

SEXP border_detection_circular_environment_cwrap(SEXP num_bins_x_r,
                                                  SEXP num_bins_y_r,
                                                  SEXP occ_map_r)
{
  int num_bins_x = INTEGER_VALUE(num_bins_x_r);
  int num_bins_y = INTEGER_VALUE(num_bins_y_r);
  double* occ_map = REAL(occ_map_r);
  
  int total_bins= num_bins_x*num_bins_y;
  
  double* tocc_map = (double*)malloc(total_bins*sizeof(double));
  for(int x = 0; x < num_bins_x; x++) 
    for(int y = 0; y < num_bins_y; y++)
      tocc_map[x*num_bins_y+y]=occ_map[x+num_bins_x*y];
  
  // 1) identify the borders in the occupancy map
  int num_bins_border=0;
  SEXP out = PROTECT(allocMatrix(INTSXP,num_bins_x,num_bins_y));
  int* border_map = (int*) malloc(total_bins*sizeof(int));
  int* border_x =  (int*) malloc(total_bins*sizeof(int)); // list of x index for border bins
  int* border_y = (int*) malloc(total_bins*sizeof(int)); // list of y index for border bins
  // Rprintf("identify_border_pixels_in_occupancy_map() call\n");
  identify_border_pixels_in_occupancy_map(tocc_map,num_bins_x,num_bins_y,border_map,border_x,border_y,&num_bins_border);
  
  
  for(int i = 0; i < num_bins_border;i++)
  {
    border_map[border_x[i]*num_bins_y+border_y[i]]=2;
  }
  
  int* o = INTEGER_POINTER(out);
  //transpose the maps for R
  for(int x =0; x < num_bins_x;x++)
    for(int y =0; y < num_bins_y;y++)
      o[y*num_bins_x+x]=border_map[x*num_bins_y+y];
  
  
  free(border_x);
  free(border_y);
  free(tocc_map);
  UNPROTECT(1);
  return (out);
  
}


SEXP border_score_circular_environment_cwrap(SEXP cells_r,
                                             SEXP cell_lines_r,
                                             SEXP num_bins_x_r,
                                             SEXP num_bins_y_r,
                                             SEXP occ_map_r,
                                             SEXP maps_r,
                                             SEXP percent_threshold_field_r,
                                             SEXP min_bins_in_field_r)
{
  int* cells = INTEGER_POINTER(cells_r);
  int cell_lines = INTEGER_VALUE(cell_lines_r);
  int num_bins_x = INTEGER_VALUE(num_bins_x_r);
  int num_bins_y = INTEGER_VALUE(num_bins_y_r);
  double* occ_map = REAL(occ_map_r);
  double* maps = REAL(maps_r);
  double* one_map;
  double* one_map_copy;
  
  int total_bins= num_bins_x*num_bins_y;
  
  // create a copy of map otherwise it will remove the fields from it
  double* maps_copy = (double*)malloc(total_bins*cell_lines*sizeof(double));
  double* tocc_map = (double*)malloc(total_bins*sizeof(double));
  // transpose the map
  for(int i =0; i < cell_lines;i++){
    one_map=maps+(i*total_bins);
    one_map_copy=maps_copy+(i*total_bins);
    for(int x = 0; x < num_bins_x; x++) 
      for(int y = 0; y < num_bins_y; y++)
        one_map_copy[x*num_bins_y+y]=one_map[x+num_bins_x*y]; 
  }
  for(int x = 0; x < num_bins_x; x++) 
    for(int y = 0; y < num_bins_y; y++)
      tocc_map[x*num_bins_y+y]=occ_map[x+num_bins_x*y];
  
  
  SEXP out = PROTECT(allocMatrix(REALSXP,5,cell_lines));
  double* o = REAL(out);
  
  double* cm = (double*) malloc(cell_lines*sizeof(double));
  double* dm = (double*) malloc(cell_lines*sizeof(double));
  double* border_score =  (double*) malloc(cell_lines*sizeof(double));
  int* num_fields_detected = (int*) malloc(cell_lines*sizeof(int));
  double* map_polarity = (double*) malloc(cell_lines*sizeof(double));
  border_score_circular_environment(cells,
                                       cell_lines,
                                       num_bins_x,
                                       num_bins_y,
                                       tocc_map,
                                       maps_copy,
                                       REAL(percent_threshold_field_r)[0],
                                       INTEGER_VALUE(min_bins_in_field_r),
                                       border_score, // border score
                                       cm,
                                       dm,
                                       num_fields_detected,
                                       map_polarity);
  
  for(int i = 0; i < cell_lines; i++)
  {
    
    if(num_fields_detected[i]>0){
      o[i*5+0]=border_score[i];
      o[i*5+1]=cm[i];
      o[i*5+2]=dm[i];
      o[i*5+3]=(double)num_fields_detected[i];
      o[i*5+4]=map_polarity[i];
    }
    else{
      o[i*5+0]=NA_REAL;
      o[i*5+1]=NA_REAL;
      o[i*5+2]=NA_REAL;
      o[i*5+3]=(double)num_fields_detected[i];
      o[i*5+4]=map_polarity[i];
    }
  }
  
  free(maps_copy);
  free(tocc_map);
  free(cm);
  free(dm);
  free(num_fields_detected);
  free(map_polarity);
  
  UNPROTECT(1);
  return (out);
  
}
void border_score_circular_environment(int* cells, 
                                       int cell_lines, 
                                       int num_bins_x, 
                                       int num_bins_y, 
                                       double* occ_map, 
                                       double* maps, 
                                       double percent_threshold_field, 
                                       int min_bins_in_field, 
                                       double* border_score, 
                                       double* cm, 
                                       double* dm, 
                                       int* num_fields_detected,
                                       double* map_polarity){
  
  // 1) identify the borders in the occupancy map
  int num_bins_border=0;
  int total_bins=num_bins_x*num_bins_y;
  int* border_map = (int*) malloc(total_bins*sizeof(int)); // map
  int* border_x =  (int*) malloc(total_bins*sizeof(int)); // list of x index for border bins
  int* border_y = (int*) malloc(total_bins*sizeof(int)); // list of y index for border bins
  identify_border_pixels_in_occupancy_map(occ_map,num_bins_x,num_bins_y,border_map,border_x,border_y,&num_bins_border);
  double* one_map;
  double max_fr;
  
  // for each valid pixel in occ_map, find the closest distance to a border pixel, get the largest closest distance
  double min_distance_one_pixel=0;
  double dist=0;
  int number_fields_detected=0;
  double* bin_distance_to_nearest_border = (double*) malloc(total_bins*sizeof(double));
  double CM=0;
  double max_CM=0;
  double DM=0;
  double WVL=0;
  int pixels_covered_one_wall;
  double* one_field_map; // rate of one firing field
  double* all_fields_map; // rate of all firing fields
  double* detection_map; // where fields get set to -1
  double max_fr_remaining=0;
  double threshold_hz;
  double mean_x_field;
  double mean_y_field;
  double radius_field;
  int num_bins_field;
  double sum_firing_rate_in_fields;
  double* rad_map;
  rad_map = (double*)malloc(total_bins*sizeof(double));
  one_field_map = (double*) malloc(total_bins*sizeof(double));
  all_fields_map = (double*) malloc(total_bins*sizeof(double));
  detection_map = (double*) malloc(total_bins*sizeof(double));
  
  double max_possible_distance_to_border=0;
  for(int x = 0; x < num_bins_x; x++)
  {
    for(int y = 0; y < num_bins_y; y++)
    {
      if (occ_map[(x*num_bins_y) + y] != -1)
      {
        min_distance_one_pixel=num_bins_x+num_bins_y;
        for (int i = 0; i < num_bins_border; i++)
        {
          dist=distance(x,y,border_x[i],border_y[i]);
          if(dist<min_distance_one_pixel)
          {
            min_distance_one_pixel=dist;
          }
        }
        if(min_distance_one_pixel>max_possible_distance_to_border)
        {
          max_possible_distance_to_border=min_distance_one_pixel;
        }
      }
    }
  }
  
 
  for(int i = 0; i < cell_lines; i++)
  { // for each cell, 
    
    max_CM=0;
    one_map = maps + (i*total_bins);	      
    for(int j = 0; j < total_bins; j++)
    {
      one_field_map[j] =-1.0;
      all_fields_map[j] =-1.0;
      detection_map[j]=one_map[j];
    }
    
    // calculate the maximum CM
    max_fr=find_max_double(total_bins,detection_map);
    max_fr_remaining=max_fr;
    threshold_hz = max_fr*percent_threshold_field/100;
    number_fields_detected=0;
    // cerr << "max fr: " << max_fr << '\n';
    while(max_fr_remaining>=threshold_hz)
    {
      detect_one_field_with_field(detection_map,
                                  num_bins_x,
                                  num_bins_y,
                                  min_bins_in_field,
                                  threshold_hz,
                                  &mean_x_field,
                                  &mean_y_field,
                                  &radius_field,
                                  &num_bins_field,
                                  -1.0,
                                  one_field_map);
      if (mean_x_field != -1)
      {// valid field detected
        number_fields_detected++;
        // copy the field to the all field maps
        for(int k=0; k < total_bins;k++)
        {
          if(one_field_map[k]!=-1.0)
            all_fields_map[k]=one_field_map[k];
        }
        
        // number of occupied wall pixels
        pixels_covered_one_wall=0;
        for(int m =0; m < num_bins_border; m++)
        {
          if(one_field_map[(border_x[m]*num_bins_y)+border_y[m]]!=-1.0)
            pixels_covered_one_wall++;
        }
        if(pixels_covered_one_wall!=0)
          CM=(double)pixels_covered_one_wall/ ((double)num_bins_border); // this is different from rectangular
        else
          CM=0;
        if(CM>max_CM)
          max_CM=CM;
      }
      max_fr_remaining=find_max_double(total_bins,detection_map);
    }
    CM=max_CM;
    
    // calculate DM
    sum_firing_rate_in_fields=sum_double(total_bins,all_fields_map,-1.0);
    DM=0;
    for(int x = 0; x < num_bins_x; x++)
    {
      for(int y = 0; y < num_bins_y; y++)
      {
        if (all_fields_map[(x*num_bins_y) + y] != -1.0)// field bin
        {
          // find the closest distance to a border
          min_distance_one_pixel=num_bins_x+num_bins_y;
          for (int j = 0; j < num_bins_border; j++)
          {
            dist=distance(x,y,border_x[j],border_y[j]);
            if(dist<min_distance_one_pixel)
            {
              min_distance_one_pixel=dist;
            }
          }
          DM+=(min_distance_one_pixel/max_possible_distance_to_border)*one_map[(x*num_bins_y) + y];
        }
      }
    }
    DM=DM/sum_firing_rate_in_fields;
    
    
    // calculate the vector length of the all_fields_map relative to the center of the occupancy map
    // weighted with the firing rate of each pixel
    
    for(int i = 0; i < total_bins;i++)
      rad_map[i]=-1.0;
    double x_center=num_bins_x/2.0;
    double y_center=num_bins_y/2.0;
    for(int x = 0; x < num_bins_x; x++)
      for(int y = 0; y < num_bins_y; y++)
        if(all_fields_map[(x*num_bins_y)+y]!=-1.0)
        {
          if(x_center-x!=0||y_center-y!=0)
            rad_map[(x*num_bins_y)+y]=direction(x_center,y_center,x,y);
        }
    //2 call a function to get a weighted vector length
    WVL=mean_vector_length_weighted(rad_map,all_fields_map,total_bins);
  
    if(number_fields_detected>0)
    {
      if(CM>1){
        Rprintf("CM is larger than 1, %lf\n",CM);
          }
      if(DM>1){
          Rprintf("DM is larger than 1, %lf\n",DM);
      }
          
      cm[i]= CM;
      dm[i]= DM;
      num_fields_detected[i]=number_fields_detected;
      border_score[i]=(CM-DM)/(CM+DM);
      map_polarity[i]=WVL;
    }
    else
    {
      cm[i]= NAN;
      dm[i]= NAN;
      num_fields_detected[i]=0;
      border_score[i]=NAN;
      map_polarity[i]=NAN;
    }
  }
  free(bin_distance_to_nearest_border);
  free(one_field_map);
  free(all_fields_map);
  free(detection_map);
  free(rad_map);
}


void spike_head_direction(double *hd_whl, 
			  int whl_lines, 
			  int *res, 
			  int res_lines, 
			  double *hd_spike, 
			  int res_samples_per_whl_sample,
			  int* start_interval, 
			  int* end_interval, 
			  int interval_lines)
{
  // get the res index for the intervals
  int* start_interval_index;
  int* end_interval_index;

  start_interval_index = (int*) malloc(interval_lines*sizeof(int));
  end_interval_index =  (int*) malloc(interval_lines*sizeof(int));
  res_index_for_intervals(&interval_lines,
			  start_interval,
			  end_interval, 
			  res_lines, 
			  res,
			  start_interval_index,
			  end_interval_index,1);
  // set all the spike position to -1
  set_array_to_value_double(hd_spike,res_lines,-1.0);

  //loop for each interval and set the position when valid//
  int whl_index;
  int residual;
  double deg_diff,corr_hd;
  for (int k = 0; k < interval_lines; k++)
    {
      for (int i = start_interval_index[k]; i < end_interval_index[k]; i++)
	{
	  
	  if (res[i] != -1.0)
	    {
	      whl_index = (res[i] / res_samples_per_whl_sample) - 1;
	      if (whl_index >=  whl_lines || (whl_index < 0 && whl_index != -1))
		{
		  Rprintf("whl file is too small for some res value\n");
		}
	      residual = res[i] % res_samples_per_whl_sample;
	      // spikes before the first frame capture
	      if (whl_index == -1)
		{hd_spike[i] = hd_whl[0];}
	      
	      // spikes between the first and last frame capture
	      if ((whl_index < whl_lines -1) && 
		  (whl_index != -1) && 
		  (hd_whl[whl_index] != -1.0) &&
		  (hd_whl[whl_index+1]!= -2))
		{
		  // circular data
		  deg_diff=sqrt((hd_whl[whl_index+1]-hd_whl[whl_index])*
				(hd_whl[whl_index+1]-hd_whl[whl_index]));
		  if (deg_diff <= 180)
		    {
		      hd_spike[i] =
			hd_whl[whl_index] +
			((hd_whl[whl_index+1]- hd_whl[whl_index])
			 /res_samples_per_whl_sample * residual);
		    }
		  else // the difference is larger than 180
		    {
		      if ( hd_whl[whl_index+1] >= hd_whl[whl_index])
			{
			  corr_hd=hd_whl[whl_index]+360;
			  hd_spike[i] =  corr_hd +
			    ((hd_whl[whl_index+1]- corr_hd)
			     /res_samples_per_whl_sample * residual);
			}
		      
		      if ( hd_whl[whl_index+1] < hd_whl[whl_index])
			{
			  corr_hd=hd_whl[whl_index+1]+360;
			  hd_spike[i] = hd_whl[whl_index] +
			    ((corr_hd- hd_whl[whl_index])
			     /res_samples_per_whl_sample * residual);
			}
		      if ( hd_spike[i] >= 360)
			{
			  hd_spike[i]= hd_spike[i]-360;
			}
		    }
		}
	      // if spikes after the last frame capture
	      if (whl_index == whl_lines -1)
		{
		  hd_spike[i] = hd_whl[whl_index];
		}
	       
	      // cope with the -1 in whl
	      if (hd_whl[whl_index] == -1.0 || hd_whl[whl_index+1]==-1.0)
		{
		  hd_spike[i] = -1.0;
		}
	    }
	  else // the res value was set at -1
	    {
	      hd_spike[i] = -1.0;
	    }
	}
    }
  free(start_interval_index);
  free(end_interval_index);
  return ;
}
SEXP spike_head_direction_cwrap(SEXP hd_whl_r,
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
  SEXP out = PROTECT(allocVector(REALSXP,res_lines));
  double* hd_spike = (double*)malloc(res_lines*sizeof(double));
  spike_head_direction(REAL(hd_whl_r),
		       whl_lines,
		       INTEGER_POINTER(res_r),
		       res_lines,
		       hd_spike,
		       INTEGER_VALUE(res_samples_per_whl_sample_r),
		       INTEGER_POINTER(start_interval_r),
		       INTEGER_POINTER(end_interval_r),
		       INTEGER_VALUE(interval_lines_r));
  double* rans;
  rans = REAL(out);
  // copy the results in a SEXP
  for(int i = 0; i < res_lines; i++) {
    rans[i] = hd_spike[i];
  }
  free(hd_spike);
  UNPROTECT(1);
  return(out);
}
void occupancy_histogram(int x_bins, // num of bins x
			 double pixels_per_bin_x,
			 double *x_whl, 
			 int whl_lines, // num lines in whl data
			 double *map, // occupancy map
			 double ms_per_sample, // ms per sample for whl data
			 int *start_interval, // starts of the intervals in res value
			 int *end_interval, // ends of intervals in res value
			 int interval_lines, // number of intervals
			 int res_samples_per_whl_sample,
			 int num_repetitions // if the occupancy size is double, then add the value
			                     // twice at value/pixels_per_bin_x and value*2/pixels_per_bin_x
			 )
{
  /**********************************************************************
     Create the occupancy histogram
     Work like occupancy_map() but for 1d data
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
  
  // set the map to 0
  for(int x = 0; x < x_bins; x++)
    {
      map[x] = 0;
    }
    
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
	  for (int j = 0; j < num_repetitions; j++)
	    {
	      bin_x = (int) x_whl[i]/pixels_per_bin_x + ((x_bins/num_repetitions)*j);
	      map[bin_x]=map[bin_x]+time_to_add;
	    }
	}
    }
  // set the bins that have less then min_occ ms
  for(int x = 0; x < x_bins; x++)
    {
      if ( map[x] <= min_occ)
	{
	  map[x] = -1.0;
	}
    }
}
SEXP occupancy_histogram_cwrap(SEXP x_bins_r, // num of bins x
			       SEXP pixels_per_bin_x_r,
			       SEXP x_whl_r, 
			       SEXP whl_lines_r, // num lines in whl data
			       SEXP ms_per_sample_r, // ms per sample for whl data
			       SEXP start_interval_r, // starts of the intervals in res value
			       SEXP end_interval_r, // ends of intervals in res value
			       SEXP interval_lines_r, // number of intervals
			       SEXP res_samples_per_whl_sample_r,
			       SEXP n_repetitions_r)
{
  int x_bins=INTEGER_VALUE(x_bins_r);
  SEXP out = PROTECT(allocVector(REALSXP,x_bins)); //  to return
  double* histo = (double*)malloc(x_bins*sizeof(double)); // to pass to c function
  occupancy_histogram(x_bins, // num of bins x
		      REAL(pixels_per_bin_x_r)[0],
		      REAL(x_whl_r),
		      INTEGER_VALUE(whl_lines_r), // num lines in whl data
		      histo, // occupancy map
		      REAL(ms_per_sample_r)[0], // ms per sample for whl data
		      INTEGER_POINTER(start_interval_r), // starts of the intervals in res value
		      INTEGER_POINTER(end_interval_r), // ends of intervals in res value
		      INTEGER_VALUE(interval_lines_r), // number of intervals
		      INTEGER_VALUE(res_samples_per_whl_sample_r),
		      INTEGER_VALUE(n_repetitions_r));


  double* rans;
  rans = REAL(out);
  // copy the results in a SEXP
  for(int x = 0; x < x_bins; x++)
    rans[x] = histo[x];
  free(histo);
  UNPROTECT(1);
  return(out);
}
void create_place_field_linear( int x_bins,
				double pixels_per_bin, 
				double *x_spike, 
				int *clu,
				int res_lines,
				int target_cell,
				double *occupancy_map,
				double *place_field,
				int num_repetitions)
{
  /*********************************************************
  function to calculate the linear place field of one cell
  *********************************************************/
  int total_bins = x_bins;
  int bin_x;
  // set the spike array to 0
  set_array_to_value_double(place_field,total_bins,0);
  
  // add the number of spikes for each bin
  for(int i = 0; i <res_lines; i++)
    {
      if (clu[i]==target_cell)
	{
	  if (x_spike[i] != -1.0)
	    {
	      for (int j = 0; j < num_repetitions; j++)
		{
		  bin_x = (int) x_spike[i]/pixels_per_bin + ((x_bins/num_repetitions)*j);
		  place_field[bin_x]++;
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
	  if (place_field[i] != 0)
	    {
	      place_field[i] = (double)place_field[i]/((double)occupancy_map[i]/1000);
	    }
	}
    }
}
SEXP firing_rate_histo_cwrap(SEXP num_bins_x_r,
			     SEXP pixels_per_bin_x_r,
			     SEXP x_spike_r, 
			     SEXP clu_r,
			     SEXP res_lines_r,
			     SEXP cells_r,
			     SEXP num_cells_r,
			     SEXP occ_histo_r,
			     SEXP smooth_histo_sd_r,
			     SEXP repetitions_r,
			     SEXP circular_smooth_r)
{

  int num_bins_x=INTEGER_VALUE(num_bins_x_r);
  int num_cells=INTEGER_VALUE(num_cells_r);
  int total_size = num_bins_x*num_cells;
  int* cells = INTEGER_POINTER(cells_r);
  int target_cell;
  SEXP out = PROTECT(allocVector(REALSXP,total_size));
  double* histos = (double*)malloc(total_size*sizeof(double)); 
  double* one_histo;
  for(int i = 0; i < num_cells;i++)
    {
      target_cell = cells[i];
      // use a pointer to place give the address for each cell
      one_histo = histos + (i*num_bins_x);
      create_place_field_linear(num_bins_x,
				                        REAL(pixels_per_bin_x_r)[0],
				                        REAL(x_spike_r),
				                        INTEGER_POINTER(clu_r),
				                        INTEGER_VALUE(res_lines_r),
				                        target_cell,
				                        REAL(occ_histo_r),
				                        one_histo,
				                        INTEGER_VALUE(repetitions_r));
      if(INTEGER_VALUE(circular_smooth_r)>0)
      {
        smooth_double_gaussian_circular(one_histo,num_bins_x, REAL(smooth_histo_sd_r)[0],-1.0);
      }
      else
      {
        smooth_double_gaussian(one_histo,num_bins_x, REAL(smooth_histo_sd_r)[0],-1.0);
      }
    }
  double* ans;
  ans=REAL(out);
  for(int i = 0; i < total_size;i++)
    ans[i]=histos[i];
  free(histos);
  UNPROTECT(1);
  return(out);
}


SEXP detect_and_remove_field_cwrap(SEXP cells_r,
                                  SEXP cell_lines_r,
                                  SEXP maps_r,
                                  SEXP x_bins_r,
                                  SEXP y_bins_r,
                                  SEXP num_fields_to_detect_r,
                                  SEXP min_num_bins_fields_r,
                                  SEXP threshold_r,
                                  SEXP invalid_r)
{
  int cell_lines =INTEGER_VALUE(cell_lines_r);
  int x_bins=INTEGER_VALUE(x_bins_r);
  int y_bins=INTEGER_VALUE(y_bins_r);
  int one_map_size=x_bins*y_bins;
  double* tmaps = (double*)malloc(one_map_size*cell_lines*sizeof(double));
  double* maps =REAL(maps_r);
  double* one_map;
  double* one_tmap;
  // transpose the maps
  for(int i = 0; i < cell_lines; i++){
    one_map=maps+(one_map_size*i);
    one_tmap=tmaps+(one_map_size*i);  
    for(int x = 0; x< x_bins; x++)
      for(int y = 0; y< y_bins; y++)
        one_tmap[x*y_bins+y]=one_map[x+x_bins*y];
  }

  for(int i = 0; i < cell_lines; i++){
    one_map=tmaps+(one_map_size*i);
    detect_and_remove_field(one_map,
                            x_bins,
                            y_bins,
                            INTEGER_VALUE(num_fields_to_detect_r),
                            INTEGER_VALUE(min_num_bins_fields_r),
                            REAL(threshold_r)[0],
                            REAL(invalid_r)[0]);

  }
  
  SEXP out = PROTECT(allocVector(REALSXP,one_map_size*cell_lines));
  maps=REAL(out);
  // transpose back to R order
  for(int i = 0; i < cell_lines;i++){
    one_tmap = tmaps + (i*one_map_size); 
    one_map = maps + (i*one_map_size); // to ship back
    for(int x =0; x < x_bins;x++)
      for(int y =0; y < y_bins;y++)
        one_map[x+x_bins*y]=one_tmap[x*y_bins+y];
  }
  
  free(tmaps);
  UNPROTECT(1);
  return(out);
}

double detect_and_remove_field(double *map,
                               int x_bins,
                               int y_bins,
                               int num_fields_to_detect,
                               int min_num_bins_fields,
                               float threshold,
                               double invalid)
{
  /********************************************************
  funtion to detect fields in a map and remove them
  **********************************************************/
  double* x_field;
  double* y_field;
  double* radius_field;
  int* num_bins_field;
  
  int index;
  x_field= (double*) malloc(num_fields_to_detect*sizeof(double));
  y_field= (double*) malloc(num_fields_to_detect*sizeof(double));
  radius_field= (double*) malloc(num_fields_to_detect*sizeof(double));
  num_bins_field= (int*) malloc(num_fields_to_detect*sizeof(int));
  
  for (int i = 0; i < num_fields_to_detect ; i++)
  {
    // that function set to invalid the fieldfield
    detect_one_field(map,
                     x_bins,
                     y_bins,
                     min_num_bins_fields,
                     threshold,
                     &x_field[i],
                     &y_field[i],
                     &radius_field[i],
                     &num_bins_field[i],
                     invalid);
  }
  for (int i = 0; i < num_fields_to_detect ; i++)
  {
    if (x_field[i] != -1 && y_field[i] != -1)
    {
      // set the mean of the field to 1 in the map
      index=get_index_from_x_and_y_bin(x_bins,y_bins,(int)x_field[i],(int)y_field[i]);
      map[index]=1;
    }
  }
  free(x_field);
  free(y_field);
  free(radius_field);
  free(num_bins_field);
  return 0;
}

SEXP autocorrelation_doughnut_cwrap(SEXP cells_r,
                                    SEXP cell_lines_r,
                                    SEXP maps_r,
                                    SEXP x_bins_r,
                                    SEXP y_bins_r,
                                    SEXP num_fields_to_detect_r,
                                    SEXP min_num_bins_fields_r,
                                    SEXP threshold_r,
                                    SEXP px_per_bin_r,
                                    SEXP invalid_r)
{
  int cell_lines =INTEGER_VALUE(cell_lines_r);
  int x_bins=INTEGER_VALUE(x_bins_r);
  int y_bins=INTEGER_VALUE(y_bins_r);
  int one_map_size=x_bins*y_bins;
  double* tmaps = (double*)malloc(one_map_size*cell_lines*sizeof(double));
  double* maps =REAL(maps_r);
  double* one_map;
  double* one_tmap;
  // transpose the maps
  for(int i = 0; i < cell_lines; i++){
    one_map=maps+(one_map_size*i);
    one_tmap=tmaps+(one_map_size*i);  
    for(int x = 0; x< x_bins; x++)
      for(int y = 0; y< y_bins; y++)
        one_tmap[x*y_bins+y]=one_map[x+x_bins*y];
  }
  
  for(int i = 0; i < cell_lines; i++){
    one_map=tmaps+(one_map_size*i);
   
   autocorrelation_doughnut(one_map,
                            x_bins,
                            y_bins,
                            REAL(px_per_bin_r)[0],
                            INTEGER_VALUE(num_fields_to_detect_r),
                            INTEGER_VALUE(min_num_bins_fields_r),
                            REAL(threshold_r)[0],
                            REAL(invalid_r)[0]);
                            
  }
  
  SEXP out = PROTECT(allocVector(REALSXP,one_map_size*cell_lines));
  maps=REAL(out);
  // transpose back to R order
  for(int i = 0; i < cell_lines;i++){
    one_tmap = tmaps + (i*one_map_size); 
    one_map = maps + (i*one_map_size); // to ship back
    for(int x =0; x < x_bins;x++)
      for(int y =0; y < y_bins;y++)
        one_map[x+x_bins*y]=one_tmap[x*y_bins+y];
  }
  
  free(tmaps);
  UNPROTECT(1);
  return(out);
}

void autocorrelation_doughnut(double* one_auto_map,
                              int auto_num_bins_x,
                              int auto_num_bins_y,
                              double pixels_per_bin,
                              int number_fields_to_detect, // min number of bins in fields
                              int min_num_bins_per_field,
                              double field_threshold, 
                              double invalid)
{
  // get the circular area that is rotated to calculate the grid score
  int num_fields = 7;
  int num_bins = auto_num_bins_x*auto_num_bins_y;
  double* original_map;
  double* x_field;
  double* y_field;
  double* radius;
  double* distance_to_mid;
  int* area_field;
  int mid_x=auto_num_bins_x/2;
  int mid_y=auto_num_bins_y/2;
  int x;
  int y;
  int d;
  double min_radius; // for the ring
  double max_radius; // for the ring
  
  original_map = (double*) malloc(sizeof(double) *num_bins);
  x_field = (double*) malloc(sizeof(double) *num_fields);
  y_field = (double*) malloc(sizeof(double) *num_fields);
  radius = (double*) malloc(sizeof(double) *num_fields);
  distance_to_mid = (double*) malloc(sizeof(double) *num_fields);
  area_field = (int*) malloc(sizeof(int) *num_fields);
  
  
  // copy the original map before detection of fields
  for (int i = 0 ; i < num_bins ; i++)
  {
    original_map[i] = one_auto_map [i];
  }
  // get the position of the 6 fields close to mid of map
  grid_closest_peaks_to_middle(one_auto_map,
                               auto_num_bins_x,
                               auto_num_bins_y,
                               pixels_per_bin,
                               number_fields_to_detect,
                               min_num_bins_per_field,
                               field_threshold,
                               invalid,
                               x_field, // res for results
                               y_field,
                               radius,
                               distance_to_mid,
                               area_field,
                               num_fields);
  
  // to make the ring that will be rotated
  min_radius=radius[0]/pixels_per_bin; 
  if (min_radius >= auto_num_bins_x/4 || min_radius >= auto_num_bins_y/4)
  {
    Rprintf("autocorrelation_doughnut\n");
    Rprintf("radius of central field was as big as the map\n");
    min_radius=min_radius/2;
  }
  max_radius=0;
  for (int i = 1 ; i < num_fields; i++)
  {  
    if (x_field[i] != -1)
    {
      if ((distance_to_mid[i]+radius[i])/pixels_per_bin > max_radius)
      {
        max_radius=(distance_to_mid[i]+radius[i])/pixels_per_bin;
      }
    }
  }
  // set the value outside the ring to invalid
  for (int i = 0 ; i < num_bins ; i++)
  {
    if (original_map[i] != invalid)
    {
      // check the distance between that bin and the center of the map
      get_x_and_y_bin_from_index(auto_num_bins_x,
                                 auto_num_bins_y,
                                 i,
                                 &x,
                                 &y);
      d=(int)distance(mid_x,mid_y,x,y);
      if (d<min_radius|| d>max_radius)
      {
        original_map[i] = invalid;
      }
    }
  }

  for(int i = 0; i < num_bins; i++)
    one_auto_map[i]=original_map[i];
  
  free(original_map);
  free(x_field);
  free(y_field);
  free(radius);
  free(distance_to_mid);
  free(area_field);
  return;
}


SEXP autocorrelation_doughnut_rotate_cwrap(SEXP cells_r,
                                    SEXP cell_lines_r,
                                    SEXP maps_r,
                                    SEXP x_bins_r,
                                    SEXP y_bins_r,
                                    SEXP num_fields_to_detect_r,
                                    SEXP min_num_bins_fields_r,
                                    SEXP threshold_r,
                                    SEXP px_per_bin_r,
                                    SEXP rotations_r,
                                    SEXP degree_r,
                                    SEXP invalid_r)
{
  int cell_lines =INTEGER_VALUE(cell_lines_r);
  int x_bins=INTEGER_VALUE(x_bins_r);
  int y_bins=INTEGER_VALUE(y_bins_r);
  int one_map_size=x_bins*y_bins;
  int rotations = INTEGER_VALUE(rotations_r);
  double* tmaps = (double*)malloc(one_map_size*cell_lines*sizeof(double));
  double* maps =REAL(maps_r);
  double* one_map;
  double* one_tmap;
  double* ptr_rotated;
  // transpose the maps
  for(int i = 0; i < cell_lines; i++){
    one_map=maps+(one_map_size*i);
    one_tmap=tmaps+(one_map_size*i);  
    for(int x = 0; x< x_bins; x++)
      for(int y = 0; y< y_bins; y++)
        one_tmap[x*y_bins+y]=one_map[x+x_bins*y];
  }
  // all cells * rotations
  double* rotated_maps = (double*) malloc(sizeof(double)*one_map_size*cell_lines*rotations);
  
  for(int i = 0; i < cell_lines; i++){
    one_map=tmaps+(one_map_size*i);
    ptr_rotated = rotated_maps+(one_map_size*rotations*i);
    autocorrelation_doughnut_rotate(one_map,
                                    x_bins,
                                    y_bins,
                                    ptr_rotated,
                                    REAL(px_per_bin_r)[0],
                                    INTEGER_VALUE(num_fields_to_detect_r),
                                    INTEGER_VALUE(min_num_bins_fields_r),
                                    REAL(threshold_r)[0],
                                    rotations,
                                    REAL(degree_r)[0],
                                    REAL(invalid_r)[0]);
    
  }
  
  SEXP out = PROTECT(allocVector(REALSXP,one_map_size*cell_lines*rotations));
  maps=REAL(out);
  // transpose back to R order
  for(int i = 0; i < cell_lines;i++)
    for(int j = 0; j < rotations; j++){
      one_tmap = rotated_maps + (i*one_map_size*rotations+ one_map_size*j); 
      one_map = maps + (i*one_map_size*rotations+ one_map_size*j) ;
      for(int x =0; x < x_bins;x++)
        for(int y =0; y < y_bins;y++)
          one_map[x+x_bins*y]=one_tmap[x*y_bins+y];
    }
  free(tmaps);
  free(rotated_maps);
  UNPROTECT(1);
  return(out);
}

void autocorrelation_doughnut_rotate(double* one_auto_map,
                              int auto_num_bins_x,
                              int auto_num_bins_y,
                              double* rotated_maps,
                              double pixels_per_bin,
                              int number_fields_to_detect, // min number of bins in fields
                              int min_num_bins_per_field,
                              double field_threshold,
                              int num_rotations,
                              double degree,
                              double invalid)
{
  // get the circular area that is rotated to calculate the grid score
  int num_fields = 7;
  int num_bins = auto_num_bins_x*auto_num_bins_y;
  double* original_map;
  double* x_field;
  double* y_field;
  double* radius;
  double* distance_to_mid;
  int* area_field;
  int mid_x=auto_num_bins_x/2;
  int mid_y=auto_num_bins_y/2;
  int x;
  int y;
  int d;
  double min_radius; // for the ring
  double max_radius; // for the ring
  double* one_map;
  original_map = (double*) malloc(sizeof(double) *num_bins);
  x_field = (double*) malloc(sizeof(double) *num_fields);
  y_field = (double*) malloc(sizeof(double) *num_fields);
  radius = (double*) malloc(sizeof(double) *num_fields);
  distance_to_mid = (double*) malloc(sizeof(double) *num_fields);
  area_field = (int*) malloc(sizeof(int) *num_fields);
  
  
  // copy the original map before detection of fields
  for (int i = 0 ; i < num_bins ; i++)
  {
    original_map[i] = one_auto_map [i];
  }
  // get the position of the 6 fields close to mid of map
  grid_closest_peaks_to_middle(one_auto_map,
                               auto_num_bins_x,
                               auto_num_bins_y,
                               pixels_per_bin,
                               number_fields_to_detect,
                               min_num_bins_per_field,
                               field_threshold,
                               invalid,
                               x_field, // res for results
                               y_field,
                               radius,
                               distance_to_mid,
                               area_field,
                               num_fields);
  
  // to make the ring that will be rotated
  min_radius=radius[0]/pixels_per_bin; 
  if (min_radius >= auto_num_bins_x/4 || min_radius >= auto_num_bins_y/4)
  {
    Rprintf("autocorrelation_doughnut\n");
    Rprintf("radius of central field was as big as the map\n");
    min_radius=min_radius/2;
  }
  max_radius=0;
  for (int i = 1 ; i < num_fields; i++)
  {  
    if (x_field[i] != -1)
    {
      if ((distance_to_mid[i]+radius[i])/pixels_per_bin > max_radius)
      {
        max_radius=(distance_to_mid[i]+radius[i])/pixels_per_bin;
      }
    }
  }
  // set the value outside the ring to invalid
  for (int i = 0 ; i < num_bins ; i++)
  {
    if (original_map[i] != invalid)
    {
      // check the distance between that bin and the center of the map
      get_x_and_y_bin_from_index(auto_num_bins_x,
                                 auto_num_bins_y,
                                 i,
                                 &x,
                                 &y);
                                 d=(int)distance(mid_x,mid_y,x,y);
                                 if (d<min_radius|| d>max_radius)
                                 {
                                   original_map[i] = invalid;
                                 }
    }
  }
  
  // create the rotated maps 
  for(int i = 0; i < num_rotations; i++)
  {
    one_map=rotated_maps+i*num_bins;
    for (int j = 0; j < num_bins; j++)
    {
      one_map[j]=original_map[j];
    }
    
    map_rotate(one_map,
               auto_num_bins_x,
               auto_num_bins_y,
               degree*i,
               invalid);
  }
  
  
  free(original_map);
  free(x_field);
  free(y_field);
  free(radius);
  free(distance_to_mid);
  free(area_field);
  return;
}
SEXP spike_position_1d_cwrap(SEXP x_whl_r,
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
      SEXP out = PROTECT(allocVector(REALSXP,res_lines));
      double* x_spike = (double*)malloc(res_lines*sizeof(double));
      spike_position_1d(REAL(x_whl_r),
                            whl_lines,
                            INTEGER_POINTER(res_r),
                            res_lines,
                            x_spike,
                            INTEGER_VALUE(res_samples_per_whl_sample_r),
                            INTEGER_POINTER(start_interval_r),
                            INTEGER_POINTER(end_interval_r),
                            INTEGER_VALUE(interval_lines_r));
      double* rans;
      rans = REAL(out);
      // copy the results in a SEXP
      for(int i = 0; i < res_lines; i++) {
        rans[i] = x_spike[i];
      }

      free(x_spike);
      UNPROTECT(1);
      return(out);
    }
    
void spike_position_1d(double *x_whl,
                    int whl_lines,
                    int *res,
                    int res_lines,
                    double *x_spike,
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
                          end_interval_index,1);
  // set all the spike position to -1
  set_array_to_value_double(x_spike,res_lines,-1.0);
  
  //loop for each interval and set the position when valid//
  int whl_index;
  int residual;
  for (int k = 0; k < interval_lines; k++)
  {
    for (int i = start_interval_index[k]; i < end_interval_index[k]; i++)
    {
      if (res[i] != -1)
      {
        whl_index = (res[i] / res_samples_per_whl_sample) - 1;
        if (whl_index >=  whl_lines || (whl_index < 0 && whl_index != -1))
        {
          Rprintf("whl file is too small for some res value\n");
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
        // cope with the -1 values from the whl file
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
  }
  free(start_interval_index);
  free(end_interval_index);
  return;
}


SEXP spike_triggered_firing_rate_maps_cwrap(SEXP num_bins_x_r,
                                            SEXP num_bins_y_r, 
                                            SEXP pixels_per_bin_x_r,
                                            SEXP pixels_per_bin_y_r,
                                            SEXP total_bins_r,
                                            SEXP cells_r,
                                            SEXP num_cells_r,
                                            SEXP x_whl_r,
                                            SEXP y_whl_r,
                                            SEXP whl_lines_r,
                                            SEXP res_r,
                                            SEXP clu_r,
                                            SEXP res_lines_r,
                                            SEXP start_interval_r,
                                            SEXP end_interval_r,
                                            SEXP interval_lines_r,
                                            SEXP ms_per_sample_r,
                                            SEXP res_samples_per_whl_sample_r,
                                            SEXP smoothing_factor_occ_r,
                                            SEXP smoothing_factor_rate_r,
                                            SEXP min_isi_ms_r,
                                            SEXP max_isi_ms_r,
                                            SEXP res_sampling_rate_r)
{
  int total_bins = INTEGER_VALUE(total_bins_r);
  int num_cells=INTEGER_VALUE(num_cells_r);
  int res_lines=INTEGER_VALUE(res_lines_r);
  int num_bins_x=INTEGER_VALUE(num_bins_x_r);
  int num_bins_y=INTEGER_VALUE(num_bins_y_r);
  
  double* x_spike = (double*)malloc(res_lines*sizeof(double));
  double* y_spike = (double*)malloc(res_lines*sizeof(double));
  double* all_occ_maps=(double*)malloc(total_bins*num_cells*sizeof(double));
  double* all_place_maps=(double*)malloc(total_bins*num_cells*sizeof(double));
  SEXP out = PROTECT(allocVector(REALSXP,total_bins*num_cells));
  
  
  create_place_firing_maps_spike_triggered(INTEGER_VALUE(num_bins_x_r),
                                           INTEGER_VALUE(num_bins_y_r),
                                           REAL(pixels_per_bin_x_r)[0], REAL(pixels_per_bin_x_r)[0],
                                           INTEGER_VALUE(total_bins_r),
                                           INTEGER_POINTER(cells_r),
                                           INTEGER_VALUE(num_cells_r),
                                           REAL(x_whl_r),
                                           REAL(y_whl_r),
                                           INTEGER_VALUE(whl_lines_r),
                                           INTEGER_POINTER(res_r),
                                           INTEGER_POINTER(clu_r),
                                           INTEGER_VALUE(res_lines_r),
                                           x_spike,
                                           y_spike,
                                           INTEGER_POINTER(start_interval_r),
                                           INTEGER_POINTER(end_interval_r),
                                           INTEGER_VALUE(interval_lines_r),
                                           all_occ_maps,
                                           all_place_maps,
                                           REAL(ms_per_sample_r)[0],
                                           INTEGER_VALUE(res_samples_per_whl_sample_r),
                                           REAL(smoothing_factor_occ_r)[0],
                                           REAL(smoothing_factor_rate_r)[0],
                                           REAL(min_isi_ms_r)[0],
                                           REAL(max_isi_ms_r)[0],
                                           INTEGER_VALUE(res_sampling_rate_r));
  double* ans;
  double* one_ans;
  double* one_map;
  ans=REAL(out);

  //transpose the maps for R
  for(int i = 0; i < num_cells;i++){
    one_map = all_place_maps + (i*total_bins);
    one_ans = ans + (i*total_bins);
    for(int x =0; x < num_bins_x;x++)
      for(int y =0; y < num_bins_y;y++)
        one_ans[x+num_bins_x*y]=one_map[x*num_bins_y+y];
  }
  free(x_spike);
  free(y_spike);
  free(all_occ_maps);
  free(all_place_maps);
  UNPROTECT(1);
  return out;
}


void create_place_firing_maps_spike_triggered(int num_bins_x,
                                              int num_bins_y, 
                                              double pixels_per_bin_x,
                                              double pixels_per_bin_y,
                                              int total_bins,
                                              int* cells,
                                              int num_cells,
                                              double* x_whl,
                                              double* y_whl,
                                              int whl_lines,
                                              int* res,
                                              int* clu,
                                              int res_lines,
                                              double* x_spike,
                                              double* y_spike,
                                              int* start_interval,
                                              int* end_interval,
                                              int interval_lines,
                                              double* all_occ_maps,
                                              double* all_place_maps,
                                              double ms_per_sample,
                                              int res_samples_per_whl_sample,
                                              double smoothing_factor_occ,
                                              double smoothing_factor_rate,
                                              double min_isi_ms,
                                              double max_isi_ms,
                                              int res_sampling_rate)
{
  int target_cell;
  double* one_place_map; // pointer
  double* one_occ_map; // pointer

  /////////////////////////////////////////////////////////////////////////////////
  ///          get the x and y position of every spikes  ///          
  ////////////////////////////////////////////////////////////////////////////////
  // outside intervals is set to -1 //
  spike_position(x_whl,
                 y_whl,
                 whl_lines,
                 res,
                 res_lines,
                 x_spike,
                 y_spike,
                 res_samples_per_whl_sample,
                 start_interval,
                 end_interval,
                 interval_lines);
  
  ///////////////////////////////////////////////////////////////////////////
  /// get spike_triggered_occupancy_map        ////
  /// unique for each cell                                   ////
  //////////////////////////////////////////////////////////////////////////
  for(int i = 0; i < num_cells; i++)
  {     
    target_cell = cells[i];
    // use a pointer to place give the address for each cell
    one_occ_map = all_occ_maps + (i*total_bins);
    spike_triggered_occupancy_map(num_bins_x,
                                  num_bins_y, 
                                  pixels_per_bin_x,
                                  pixels_per_bin_y,
                                  x_whl,
                                  y_whl,
                                  whl_lines,
                                  one_occ_map,
                                  ms_per_sample,
                                  start_interval,
                                  end_interval,
                                  interval_lines,
                                  res_samples_per_whl_sample,
                                  res,
                                  clu,
                                  res_lines,
                                  target_cell,
                                  min_isi_ms,
                                  max_isi_ms,
                                  res_sampling_rate,
                                  x_spike,
                                  y_spike);
    
    // smooth the occupancy map
    smooth_double_gaussian_2d(one_occ_map,num_bins_x, num_bins_y, smoothing_factor_occ,-1.0);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////
  ///       get the spike triggered place field for every cell    ////
  //////////////////////////////////////////////////////////////////////////////////////////////
  
  for(int i = 0; i < num_cells; i++)
  {     
    target_cell = cells[i];
    // use a pointer to place give the address for each cell
    one_place_map = all_place_maps + (i*total_bins);
    one_occ_map = all_occ_maps + (i*total_bins);
    
    spike_triggered_place_map(num_bins_x,
                              num_bins_y, 
                              pixels_per_bin_x,
                              pixels_per_bin_y,
                              x_spike, 
                              y_spike,
                              res,
                              clu,
                              res_lines,
                              target_cell,
                              target_cell,
                              one_occ_map,
                              one_place_map,
                              min_isi_ms,
                              max_isi_ms,
                              res_sampling_rate);
    
    /// smooth the place map
    smooth_double_gaussian_2d(one_place_map,num_bins_x, num_bins_y, smoothing_factor_rate,-1.0);
  }
}



void spike_triggered_occupancy_map(int x_bins, // num of bins x
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
                                   int res_samples_per_whl_sample,
                                   int* res,
                                   int* clu,
                                   int res_lines,
                                   int target_cell,
                                   double min_isi_ms,
                                   double max_isi_ms,
                                   int res_sampling_rate,
                                   double* x_spike,
                                   double* y_spike)
{
  /*****************************************************************************
  Create the occupancy map for 2d firing map with spike triggered method
  
  This add  ms_per_sample for every valid whl value that is in 
  the intervals given valid whl values outside the intervals 
  will not be added to occupancy map.
  If a valid whl value is for a period outside intervals,
  the right proportion of time will be removed
  
  There is a minimum time that the animal has to 
  spend in a bin for it to be valid; see variable min_occ
  
  Bins that are on there own in space are removed.
  ********************************************************************************/
  int min_occ = 500; // minimum time spent in a bin
  // variable to deal with the intervals
  int res_value_for_whl_sample; // time in res associated with a whl sample
  int start_res_value; // starting time of a whl sample in res value
  int end_res_value; // end time of a whl sample
  int res_value_to_add; // number of res inside the interval
  int interval_index=0;
  double time_to_add;
  int bin_x; // where to add the time in the map
  int bin_y;
  int res_start_window; // start of window around the spike
  int res_end_window; // end of window around the spike
  int whl_index_start_window; // index in whl array of the start of window around the spike
  int whl_index_end_window; // index in whl array of the end of window around the spike
  double relative_x_whl;
  double relative_y_whl;
  int min_isi_res=(int)(min_isi_ms*res_sampling_rate/1000);
  int max_isi_res=(int)(max_isi_ms*res_sampling_rate/1000);
  
  
  // set the map to 0
  for(int x = 0; x < x_bins; x++)
  {
    for(int y = 0; y < y_bins; y++)
    {
      map[(x*y_bins) + y] = 0;
    }
  }
  
  // loop for each valid spike of the target cell (x_spike and y_spike != -1)
  for (int i = 0; i < res_lines; i++)
  {
    if(clu[i]==target_cell&&res[i]!=-1&&x_spike[i]!=-1&&y_spike[i]!=-1)
    {
      // find start and end of window
      res_start_window=res[i]+min_isi_res;
      res_end_window=res[i]+max_isi_res;
      
      if(res_start_window<0) res_start_window=0; // make sure this is not set to negative value
      
      // get index of whl for the start and end of window
      whl_index_start_window=(res_start_window/res_samples_per_whl_sample)-1;
      whl_index_end_window=(res_end_window/res_samples_per_whl_sample)-1;
      if (whl_index_start_window>= whl_lines) whl_index_start_window=whl_lines-1;
      if (whl_index_start_window<0) whl_index_start_window=0;
      if (whl_index_end_window>=whl_lines) whl_index_end_window=whl_lines-1;
      if (whl_index_end_window<0) whl_index_end_window=0;
      
      // add the time of the window around this spike
      for(int j = whl_index_start_window; j < whl_index_end_window; j++)
      {
        if (x_whl[j] != -1.0)
        {
          // get the begining of the period cover by this whl sample, in res value
          res_value_for_whl_sample=(res_samples_per_whl_sample*j)+res_samples_per_whl_sample;
          start_res_value=res_value_for_whl_sample-(res_samples_per_whl_sample/2); // start of whl sample
          end_res_value=res_value_for_whl_sample+(res_samples_per_whl_sample/2); // end of whl sample
          
          // start with the assumption that time to add is  = 0; as if outside the valid intervals
          res_value_to_add=0;
          
          // then add res_time if there are intervals around start_res_value and end_res_value
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
          
          if(res_value_to_add>0)
          {
            // calculate the time in ms, and add it to the right bin of the occ map
            time_to_add=ms_per_sample*((double)res_value_to_add/res_samples_per_whl_sample);
            
            // get the relative position of whl for this sample
            relative_x_whl=x_whl[j]-x_spike[i];
            relative_y_whl=y_whl[j]-y_spike[i];
            
            
            bin_x = (int)((x_bins/2)+(relative_x_whl/pixels_per_bin_x));
            bin_y = (int)((y_bins/2)+(relative_y_whl/pixels_per_bin_y));
            
            // check if it falls in the map
            if ((bin_x<x_bins)&&(bin_y<y_bins))
            {
              map[(bin_x * y_bins) + bin_y]=
                map[(bin_x * y_bins) + bin_y] +
                time_to_add;
            }
          }
        }
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
void spike_triggered_place_map(int x_bins,
                               int y_bins, 
                               double pixels_per_bin_x, 
                               double pixels_per_bin_y, 
                               double *x_spike, 
                               double *y_spike, 
                               int* res, 
                               int *clu, 
                               int res_lines, 
                               int target_cell1, 
                               int target_cell2,
                               double *occupancy_map, 
                               double *place_map, 
                               double min_isi_ms,
                               double max_isi_ms,
                               int res_sampling_rate)
{
  /*******************************************************
  function to calculate the firing rate map of one cell
  needs an occupancy map and the x and y position
  of all the spikes
  
  the map is constructed relative to each spike.
  *********************************************************/
  int total_bins = x_bins * y_bins;
  int bin_x;
  int bin_y;
  int j;
  double relative_x, relative_y;
  int res_start_window, res_end_window;
  int min_isi_res=(int)(min_isi_ms*res_sampling_rate/1000);
  int max_isi_res=(int)(max_isi_ms*res_sampling_rate/1000);
  
  
  // set the place_field array to 0
  for (int i = 0; i < total_bins ; i++)
  {
    place_map[i] = 0;
  }
  
  // add the number of spikes in the place_field array
  for(int i = 0; i < res_lines; i++)
  {
    if (clu[i]==target_cell1&& res[i]!=-1&&x_spike[i]!=-1.0)
    {
      
      // add the spikes around this one that are within the time period
      res_start_window=res[i]+min_isi_res;
      res_end_window=res[i]+max_isi_res;
      
      
      j=i+1;
      while(res[j]<res_end_window&&j<res_lines)
      {
        // if spike of the same cell and res and position are valid
        if(clu[j]==target_cell2&&res[j]>res_start_window&&res[j]!=-1.0&&x_spike[j]!=-1.0)
        {
          // get the relative position of the spike
          relative_x=x_spike[j]-x_spike[i];
          relative_y=y_spike[j]-y_spike[i];
          // add spike, 0,0 is the middle of the map
          bin_x = (int)((x_bins/2)+(relative_x/pixels_per_bin_x));
          bin_y = (int)((y_bins/2)+(relative_y/pixels_per_bin_y));
          
          // check if it falls in the map
          if ((bin_x<x_bins)&&(bin_y<y_bins))
          {
            place_map[(bin_x * y_bins) + bin_y]++;
          }
        }
        j++;
      }
    }
  }
  ////////////// get the firing rate for every bin////////////////
  for (int i = 0; i < total_bins; i++)
  {
    if (occupancy_map[i] == -1.0)
    {
      place_map[i] = -1.0;
    }
    
    if (occupancy_map[i] != -1.0)
    {
      if (place_map[i] != 0) // 
      {
        place_map[i] = (double)place_map[i]/((double)occupancy_map[i]/1000);
      }
    }
  }
}



SEXP maps_rotate_cwrap(SEXP maps_r,
                       SEXP num_bins_x_r,
                       SEXP num_bins_y_r, 
                       SEXP num_cells_r,
                       SEXP degree_r)
{
  
  // transpose the maps
  double* maps = REAL(maps_r);
  int num_cells= INTEGER_VALUE(num_cells_r);
  double degrees = REAL(degree_r)[0];
  int num_bins_x = INTEGER_VALUE(num_bins_x_r);
  int num_bins_y = INTEGER_VALUE(num_bins_y_r);
  int total_bins= num_bins_x*num_bins_y;
  
  
  // create a copy with transposed maps
  double* all_maps_copy = (double*)malloc(total_bins*num_cells*sizeof(double));
  double* one_map;
  double* one_map_copy;
  
  for(int i = 0; i < num_cells; i++){
    one_map = maps + (i*total_bins);
    one_map_copy = all_maps_copy + (i*total_bins);
    for(int x =0; x < num_bins_x;x++)
      for(int y =0; y < num_bins_y;y++)
        one_map_copy[x*num_bins_y+y]=one_map[x+num_bins_x*y]; // needs to be transpose for c code
  }
  // rotate the maps
  for(int i = 0; i < num_cells; i++){
    one_map_copy = all_maps_copy + (i*total_bins);
    map_rotate(one_map_copy,
              num_bins_x,
              num_bins_y,
              degrees,
              -1.0);
  }
  SEXP out = PROTECT(allocVector(REALSXP,total_bins*num_cells));
  double* ans;
  double* one_ans;
  ans=REAL(out);

  //transpose the maps for R
  for(int i = 0; i < num_cells;i++){
    one_map = all_maps_copy + (i*total_bins);
    one_ans = ans + (i*total_bins);
    for(int x =0; x < num_bins_x;x++)
      for(int y =0; y < num_bins_y;y++)
        one_ans[x+num_bins_x*y]=one_map[x*num_bins_y+y];
  }
  
  // return the maps
  UNPROTECT(1);
  return out;
}