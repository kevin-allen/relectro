#include <R.h>
#include <Rdefines.h>

void resIndexForIntervals(int* interval_lines, // int by reference
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
  int max_start=0;
  int max_end=0;
  int invalid_int;
  int pre_end_index;
  int start_looping_index;

  for(int i = 0; i < *interval_lines;i++)
    if(start[i]>max_start)
      max_start=start[i];

  for(int i = 0; i < *interval_lines;i++)
    if(end[i]>max_end){
      max_end=end[i];
    }

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
	      printf("interval after the max res value\n");
	      invalid_int++;
	    }
	}
      (*interval_lines)=*interval_lines-invalid_int;
      if (invalid_int>0)
	{
	  printf("%d intervals have been removed",invalid_int);
	}
    }

  // that is chronology between and within interval is assumed, speed up the loop
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
 

  resIndexForIntervals(&interval_lines, // by reference, address
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
    rans[i] = start_index[i];
    rans[i+interval_lines]=end_index[i];
  }

  //free memory for rate
  free(start_index);
  free(end_index);
  UNPROTECT(4);
  return(ans);
}

