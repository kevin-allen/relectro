#include "relectro.h"

int check_interval_chronology_between(int num_lines, 
				      int* start, 
				      int* end)
{ 
  /* returns 1 if the intervals are in chronological order
     only between is checked
  */
  
  if(num_lines==1)
    return 1;
  if(num_lines==0)
    return 0;
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
			     int* end_interval_index,
			     int remove_interval_after_last_res)
{
  /* gives the res index of the first and last spikes within the intervals*/
  // the values at start_interval_index and end_interval_index are within the intervals
  // If no spike are withing the intervals, the indices are set to out of bound values
  // intervals need to be chronologically organized
  int max_res=res[res_lines-1];
  int max_start=find_max(*interval_lines, start);
  int max_end=find_max(*interval_lines, end);
  int invalid_int;
  int pre_end_index;
  int start_looping_index;
  int valid_found;
  
  
  
  
  if(remove_interval_after_last_res==1){
  // check that the intervals are in the recorded time
  if ((max_start > max_res + 1)||(max_end>max_res+1))
  {
    invalid_int=0;
    for (int i = 0; i < *interval_lines;i++)
    {
      if ((start[i]< max_res)&&(end[i]>max_res))
      {
        // this was recently comment out to have end interval later than last spike when needed.
        //end[i]=max_res+1; // will include the last spike in interval
      }
      if (start[i]> max_res) // whole interval after last spike, remove as the indices would not make sense
      {
        invalid_int++;
      }
    }
    (*interval_lines)=(*interval_lines)-invalid_int;
    if (invalid_int>0)
    {
      // this assumes that the intervals are chronologically ordered
      Rprintf("start of an interval is larger than largest res value\n");
      Rprintf("%d intervals have been eliminated\n",invalid_int);
    }
  }
  }
  /// depending on whether the intervals are chronologically organized /// 
  if(check_interval_chronology_between(*interval_lines, 
                                       start, 
                                       end)==0)
  {
    //Rprintf("no chrono\n");
    // chronology between interval not assumed, slows things down a bit
    for(int i = 0; i < *interval_lines; i++)
    { 
      start_interval_index[i]=0;
      end_interval_index[i]=res_lines-1; // by default, the last spike is valid
      valid_found=0;
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
        if(res[j] > start[i]) // needs to be within, not on interval limit
        {
          start_interval_index[i] = j;
          j = res_lines;
          valid_found++;
        }
      }
      for(int j=start_interval_index[i]; j < res_lines; j++)
      {
        if(res[j]>= end[i]) // needs to be within, not on interval limit
        {
          end_interval_index[i] = j-1;
          j = res_lines;
        }
      }
      if(valid_found==0){ // set to invalid values
        start_interval_index[i]=res_lines;
        end_interval_index[i]=res_lines;
      }
      
    }
  }
  else // that is chronology between and within interval is assumed, speed up the loop
  {
  //  Rprintf("chrono\n");
    for(int i = 0; i < *interval_lines; i++)
    { 
      start_interval_index[i]=0;
      end_interval_index[i]=res_lines-1;
      valid_found=0;
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
        if(res[j] > start[i]) // needs to be within, not on interval limit
        {
          start_interval_index[i] = j;
          j = res_lines;
          valid_found++;
        }
      }
      // find the end index
      for(int j= pre_end_index; j < res_lines; j++)
      {
        if(res[j]>=end[i]) // needs to be within, not on interval limit
        {
          end_interval_index[i] = j-1;
          j = res_lines;
        }
      }
      
      if(valid_found==0){ // set to invalid values
        start_interval_index[i]=res_lines;
        end_interval_index[i]=res_lines;
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
				SEXP res_r,
				SEXP remove_intervals_after_last_res_r)

{ 
  // transform the SEXP object in their correct c types
  // create the list of variable of correct c types


  int interval_lines;
  int* start; 
  int* end; 
  int res_lines; 
  int* res; 
  int remove_intervals_after_last_res;
  
  /* // coersion */
  interval_lines=INTEGER_VALUE(interval_lines_r);
  start=INTEGER_POINTER(start_r);
  end=INTEGER_POINTER(end_r);
  res=INTEGER_POINTER(res_r);
  res_lines=INTEGER_VALUE(res_lines_r);
  remove_intervals_after_last_res=INTEGER_VALUE(remove_intervals_after_last_res_r);
  
  if(remove_intervals_after_last_res!=0&&remove_intervals_after_last_res!=1){
    Rprintf("remove_intervals_after_last_res is not 0 or 1, set to 0\n");
    remove_intervals_after_last_res=0;
  }
  
  
  
  int* start_index = (int*)malloc(interval_lines*sizeof(int));
  int* end_index = (int*)malloc(interval_lines*sizeof(int));
 

  res_index_for_intervals(&interval_lines, // by reference, address
			  start, 
			  end, 
			  res_lines, 
			  res, 
			  start_index,
			  end_index,
			  remove_intervals_after_last_res);

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
  UNPROTECT(1);
  return(ans);
}

SEXP resWithinIntervals(SEXP interval_lines_r,
			SEXP start_r, 
			SEXP end_r, 
			SEXP res_lines_r, 
			SEXP res_r)
{

  int interval_lines = INTEGER_VALUE(interval_lines_r); ;
  int res_lines = INTEGER_VALUE(res_lines_r);
  int* start = INTEGER_POINTER(start_r);
  int* end = INTEGER_POINTER(end_r);
  int* res = INTEGER_POINTER(res_r);
  SEXP results;

  //  Rprintf("res_lines: %d, interval_lines: %d\n",res_lines,interval_lines);
  PROTECT(results = allocVector(INTSXP,res_lines));
  
  for(int i = 0; i < res_lines; i++){
    INTEGER(results)[i]=0; // out by default
  }

  for(int i = 0; i < res_lines; i++)
    for(int j = 0; j < interval_lines;j++){
      //    Rprintf("res index: %d, interval index: %d\n",i,j);
      if(res[i]>=start[j]&&res[i]<=end[j])
	{
	  INTEGER(results)[i]=1;
	  j=interval_lines;
	}
    }
  UNPROTECT(1);
  return(results);
}



void join_adjacent_intervals(int* start,
                             int* end,
                             int* num) // num of intervals
{
  /* joins adjacent intervals (end = next start) into single intervals
  warning : the arrays and num will be modified
  */

  int index=0;
  int same, j, k;
  for (int i = 0; i < *num; i++)
  {
    same=1;
    j=i+1;
    k=i;
    // check if next interval has the same start as 
    while((same==1) && (j < *num))
    {
      if (end[k] == start[j])
      {
        same=1;
        j++;
        k++;
      }
      else
      {
        same=0;
      }
    }
    start[index] = start[i];
    end[index] =  end[j-1];
    i=j-1;
    index++;
  }
  *num=index; // to return the number of intervals after joining
  return;
}

int num_intervals_after_joining_AND(int* start_1,
                                    int* end_1,
                                    int num_1,
                                    int* start_2,
                                    int* end_2,
                                    int num_2)
{
  /* gives the number of intervals that would result
  if two lists of intervals were joined in a AND fashon
  this means that a period needs to be covered by both series of intervals
  to be there at the end
  So the number is the number of time two intervals overlap
  */
 
  int num_intervals = 0;
  for (int i = 0; i < num_1; i++)
  {
    for (int j = 0 ; j < num_2; j++)
    {
      // check if there is overlap between the two intervals
      if ((start_2[j]>=start_1[i] && start_2[j] < end_1[i]) ||
          (start_1[i]>=start_2[j] && start_1[i] < end_2[j]))
      {
     //   Rprintf("i: %d, j: %d, start_1: %d, end_1: %d, start_2: %d, end_2: %d\n",i,j,start_1[i],end_1[i],start_2[j],end_2[j]);
        num_intervals++;
      }
    }
  }
  return num_intervals;
}
void join_two_lists_of_intervals_AND(int* start_1,
                                     int* end_1,
                                     int num_1,
                                     int* start_2,
                                     int* end_2,
                                     int num_2,
                                     int* start_res,
                                     int* end_res)
{
  /* joins intervals using and AND logic
  this means that a period needs to be covered by both series of intervals
  to be there at the end
  */
  
  int index=0;
  for (int i = 0; i < num_1; i++)
  {
    for (int j = 0 ; j < num_2; j++)
    {
      // check if there is overlap between the two intervals
      if ((start_2[j]>=start_1[i] && start_2[j] < end_1[i]) ||
          (start_1[i]>=start_2[j] && start_1[i] < end_2[j]))
      {
        // find the beginning of overlap, largest start
        if (start_1[i]>= start_2[j])
        {
          start_res[index]=start_1[i];
        }
        else
        {
          start_res[index]=start_2[j];
        }
        // find the end of overlap, smallest end
        if (end_1[i]<= end_2[j])
        {
          end_res[index]=end_1[i];
        }
        else
        {
          end_res[index]=end_2[j];
        }
        index++;
      }
    }
  }
  return;
}



SEXP joinIntervalsAND_cwrap(SEXP s1_r,
                            SEXP e1_r, 
                            SEXP l1_r,
                            SEXP s2_r,
                            SEXP e2_r,
                            SEXP l2_r)
{
  int* s1=INTEGER_POINTER(s1_r);
  int* e1=INTEGER_POINTER(e1_r);
  int l1=INTEGER_VALUE(l1_r);
  int* s2=INTEGER_POINTER(s2_r);
  int* e2=INTEGER_POINTER(e2_r);
  int l2=INTEGER_VALUE(l2_r);
  int* s3;
  int* e3;
  int l3;
  
  join_adjacent_intervals(s1, e1, &l1);
  join_adjacent_intervals(s2, e2, &l2);
  //Rprintf("Number of intervals after join_adjacent_intervals, l1: %d l2: %d\n",l1,l2);
  l3=num_intervals_after_joining_AND(s1, e1,l1, s2, e2, l2);
  //Rprintf("num_intervals_after_joining_AND: %d\n",l3);
  if(l3<1) return(R_NilValue);
  s3 = (int*) malloc(sizeof(int)*l3);
  e3 = (int*) malloc(sizeof(int)*l3);
  join_two_lists_of_intervals_AND(s1, e1, l1,
                                  s2, e2, l2,
                                  s3, e3);
  SEXP out = PROTECT(allocMatrix(INTSXP,l3,2));
  int* ptr=INTEGER_POINTER(out);
  for(int i = 0; i < l3;i++){
    ptr[i]=s3[i];
    ptr[i+l3]=e3[i];
  }
  free(s3);
  free(e3);
  UNPROTECT(1);
  return(out);
}
