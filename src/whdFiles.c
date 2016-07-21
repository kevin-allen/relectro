#include "relectro.h"
SEXP whd_file(SEXP x_r,
			SEXP y_r, 
			SEXP hd_r, 
			SEXP up_r,
			SEXP len_r,
			SEXP max_res_r, 
			SEXP res_per_whd_r,
			SEXP max_up_diff_r)
{
  int max_res = INTEGER_VALUE(max_res_r); ;
  int res_per_whd = INTEGER_VALUE(res_per_whd_r);
  int max_up_diff = INTEGER_VALUE(max_up_diff_r);
  double* x = REAL(x_r);
  double* y = REAL(y_r);
  double* hd = REAL(hd_r);
  double len = INTEGER_VALUE(len_r);
  
  int* up = INTEGER_POINTER(up_r);
  
  // length of whd
  int lwhd=(max_res/res_per_whd);
  
  SEXP out;
  PROTECT(out = allocMatrix(REALSXP,lwhd,3));
  double* ptr;
  ptr = REAL(out);
  
  int index_larger = 0;
  int index_smaller = 0;
  int previous_index = 0;
  double position_x,position_y,deg;
  double time_diff;
  double position_diff;
  double proportion;
  double hd_diff;
  int whd_index=0;
  
  
  for (int i = res_per_whd; i < max_res; i=i+res_per_whd)
  {
    // find the index of the up for the first value larger than i
    index_larger=0;
    index_smaller=0;
    position_x=-1;
    position_y=-1;
    deg=-1;
    
    for (int j = previous_index; j < len ; j++)
    {
      if (up[j]>i)
      {
        index_larger=j;
        index_smaller=j-1;
        previous_index=j; // to save time on next sample
        j = len; // to exit inner for loop
      }
    }
    
    // if the index_larger is 0 then -1  -1
    // that means that before tracking starts or that it is after it ends
    if (index_larger == 0)
    {
      position_x=-1;
      position_y=-1;
      deg=-1;
    }
    // if the ups were really far apart for some reason, fill with -1 (when creating the main whl file)
    else if(up[index_larger]-up[index_smaller]>max_up_diff)
    {
      position_x=-1;
      position_y=-1;
      deg=-1;
    }
    else
    {	
      // if one invalid, set to invalid
      if (x[index_smaller] != -1.0 && x[index_larger] != -1.0 && hd[index_smaller] !=-1.0 && hd[index_larger] != -1.0)
      {
        time_diff= up[index_larger] - up[index_smaller];
        proportion=(i-up[index_smaller])/time_diff;
        position_diff= x[index_larger] - x[index_smaller];
        position_x= x[index_smaller] + position_diff * proportion;
        position_diff= y[index_larger] - y[index_smaller];
        position_y= y[index_smaller] + position_diff * proportion;
        
        hd_diff=hd[index_larger]-hd[index_smaller];
        if (hd_diff<-180)
        {hd_diff=(hd[index_larger]+360)-hd[index_smaller];}
        if (hd_diff>180)
        {hd_diff=hd[index_larger]-(hd[index_smaller]+360);}
        deg=hd[index_smaller] + hd_diff * proportion;
        if (deg<0)
        { deg=deg+360;}
        if (deg>360)
        {deg=deg-360;}
      }
      else 
      {
        position_x=-1;
        position_y=-1;
        deg=-1;
      }
    } 
    
    ptr[whd_index]=position_x;
    ptr[lwhd*1+whd_index] = position_y;
    ptr[lwhd*2+whd_index] = deg;
    whd_index++;
  }
  
  
  UNPROTECT(1);
  return(out);
}


