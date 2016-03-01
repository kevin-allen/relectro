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

  spike_position(REAL(speed_r),
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

void spike_position(double *x_whl,
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
