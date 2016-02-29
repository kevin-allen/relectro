#include "relectro.h"

int check_interval_chronology_between(int num_lines, 
				      int* start, 
				      int* end)
{ 
  /* returns 1 if the intervals are in chronological order
     only between is checked
  */
  for(int i = 1; i < num_lines; i++)
    {
      if(start[i] < end[i-1])
	{      
	  return 0;
	}
    }
  return 1;
}

void res_index_for_intervals(int* interval_lines,
			     int* start, 
			     int* end, 
			     int res_lines, 
			     int* res, 
			     int* start_interval_index,
			     int* end_interval_index)
{
  /* gives the res index for the beginning and end of each intervals*/
  // the values at start_interval_index and end_interval_index are within the intervals
  // if not chronologically sorted, that takes more time
  int max_res=res[res_lines-1];
  int max_start=find_max(*interval_lines, start);
  int max_end=find_max(*interval_lines, end);
  int invalid_int;
  int pre_end_index;
  int start_looping_index;

  // check that the intervals are in the recorded time
  if ((max_start > max_res + 1)||(max_end>max_res+1))
    {
      invalid_int=0;
      for (int i = 0; i < *interval_lines;i++)
	{
	  if ((start[i]< max_res)&&(end[i]>max_res))
	    {
	      end[i]=max_res;
	    }
	  if (start[i]> max_res)
	    {
	
	      printf("invalid interval because start is larger than largest res value\n");
	      invalid_int++;
	    }
	}
      (*interval_lines)=(*interval_lines)-invalid_int;
      if (invalid_int>0)
	{
	  // this assumes that the intervals are chronologically ordered
	  printf("start of an interval is larger than largest res value\n");
	  printf("%d intervals have been eliminated\n",invalid_int);
	}
    }


  /// depending on whether the intervals are chronologically organized /// 
  if(check_interval_chronology_between(*interval_lines, 
				       start, 
				       end)==0)
    {
      // chronology between interval not assumed, slows things down a bit
      for(int i = 0; i < *interval_lines; i++)
	{ 
	  start_interval_index[i]=0;
	  end_interval_index[i]=res_lines-1; // should this be res_lines-1, was res_lines, which is out of bound of the res array?
	  
	  // find at which approximate res index you should start to loop
	  if (start[i]<res[res_lines/2])
	    {
	      if (start[i]<res[res_lines/4])
		{
		  start_looping_index=0;
		}
	      else
		{
		  start_looping_index=res_lines/4;
		}
	    }
	  else
	    {
	      if (start[i]>=res[(res_lines*3)/4])
		{
		  start_looping_index=(res_lines*3)/4;
		}
	      else
		{
		  start_looping_index=res_lines/2;
		}
	    }
	  
	  for(int j=start_looping_index; j < res_lines; j++)
	    {
	      if(res[j] >= start[i])
		{
		  start_interval_index[i] = j;
		  j = res_lines;
		}
	    }
	  for(int j=start_interval_index[i]; j < res_lines; j++)
	    {
	      if(res[j]>= end[i])
		{
		  end_interval_index[i] = j-1;
		  j = res_lines;
		}
	    }
	}
    }
  else // that is chronology between and within interval is assumed, speed up the loop
    {
      for(int i = 0; i < *interval_lines; i++)
	{ 
	  start_interval_index[i]=0;
	  end_interval_index[i]=res_lines;
	  if (i==0)
	    {	 
	      pre_end_index=0;
	    }
	  else
	    {
	      pre_end_index=end_interval_index[i-1];
	    }
	  // find the start index
	  for(int j= pre_end_index; j < res_lines; j++)
	    {
	      if(res[j] > start[i])
		{
		  start_interval_index[i] = j;
		  j = res_lines;
		}
	    }
	  // find the end index
	  for(int j= pre_end_index; j < res_lines; j++)
	    {
	      if(res[j]>=end[i])
		{
		  end_interval_index[i] = j-1;
		  j = res_lines;
		}
	    }
	}
    }
  return;
}








// wrapper to call resIndexForIntervals from .Call()
SEXP resIndexForIntervals_cwrap(SEXP interval_lines_r,
				SEXP start_r, 
				SEXP end_r, 
				SEXP res_lines_r, 
				SEXP res_r)

{ 
  // transform the SEXP object in their correct c types
  // create the list of variable of correct c types


  int interval_lines;
  int* start; 
  int* end; 
  int res_lines; 
  int* res; 

  
  // protect R object created in c code so that R does not delete them
  PROTECT(start_r=AS_INTEGER(start_r));
  PROTECT(end_r=AS_INTEGER(end_r));
  PROTECT(res_r=AS_INTEGER(res_r));
  
  /* // coersion */
  interval_lines=INTEGER_VALUE(interval_lines_r);
  start=INTEGER_POINTER(start_r);
  end=INTEGER_POINTER(end_r);
  res=INTEGER_POINTER(res_r);
  res_lines=INTEGER_VALUE(res_lines_r);
  
  
  int* start_index = (int*)malloc(interval_lines*sizeof(int));
  int* end_index = (int*)malloc(interval_lines*sizeof(int));
 

  res_index_for_intervals(&interval_lines, // by reference, address
			  start, 
			  end, 
			  res_lines, 
			  res, 
			  start_index,
			  end_index);

  // prepare a SEXP matrix to store the results
  SEXP ans;  //                       /nx /ny
  double* rans;
  PROTECT(ans = allocMatrix(REALSXP, 2, interval_lines));
  rans = REAL(ans);

  // copy the results in a SEXP
  for(int i = 0; i < interval_lines; i++) {
    rans[i*2] = start_index[i];
    rans[i*2+1]=end_index[i];
  }

  //free memory for rate
  free(start_index);
  free(end_index);
  UNPROTECT(4);
  return(ans);
}

