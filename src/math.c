#include <math.h>
#include "relectro.h"

int find_max(int num_data, int* data)
{
  /* returns the maximum value in an array */
  int max=data[0];
  for (int i =1; i < num_data; i++)
    {
      if (data[i]>max)
	{
	  max=data[i];
	}
    }
  return max;
}
void gaussian_kernel(double* kernel,
		     int size,
		     double standard_deviation)
{
  /* function to build a 1d gaussian kernel */
  if (size%2==0)
    {
      printf("size of gaussian kernel should be an odd number\n");
      return;
    }
  if (standard_deviation<=0)
    {
      printf("standard deviation for gaussian_kernel is <= 0\n");
      return;
    }
  int x;
  double part_1;
  double part_2;
  double num_part_2;
  double den_part_2;
  part_1=1.0/sqrt(2*M_PI*pow(standard_deviation,2)); // part 1 of equation do not change
  
  for (int i = 0; i < size; i++)
    {
      x=i-(size-1)/2; // distance from middle point
      num_part_2=pow(x,2);
      den_part_2=2*pow(standard_deviation,2);
      part_2= exp(-(num_part_2/den_part_2));
      kernel[i]=part_1*part_2;
    }
}
void smooth_double_gaussian(double* array, int array_size, double smooth, double invalid)
{
  /* Smooth a double array of size "arraysize" using a Gaussian kernal
     Arguments:
        double *array         pointer to data to be smoothed (memory must be pre-allocated)
        int size              number of elements in original array
        int smooth            standard deviation of the kernel
        int invalid           invalid data value to be excluded from smoothing

	the gaussian kernel has a size of 3 standard deviation on both side of the middle
  */
  if(array_size<=0)
    {
      printf("in smooth_double_gaussian size of data <=0\n");
      return;
    }
  if(smooth==0)
    {
      return; // do nothing 
    }
  if(smooth<=0)
    {
      printf("in smooth_double_gaussian standard deviation of gaussian kernel is <= 0\n");
      return;
    }
  int num_standard_deviations_in_kernel=3;
  int kernel_size=((int)(smooth*num_standard_deviations_in_kernel*2))+1; // for both side
  if (kernel_size%2!=1) // should be an odd number
    {
      kernel_size++;
    }
  double* kernel;
  double* results;
  double sum_weight;
  double sum_value;
  int index_value_for_kernel;
  // make a gaussian kernel

  kernel= (double*)malloc(kernel_size*sizeof(double));
  results = (double*)malloc(array_size*sizeof(double));
  gaussian_kernel(kernel,
		  kernel_size,
		  smooth);// standard deviation in kernel
  // for each bin of the array
  for(int i = 0; i < array_size; i++)
    {
      sum_weight=0;
      sum_value=0;
      // do the sum of kernel_weight*value_array
      for (int j = 0; j < kernel_size; j++)
	{
	  index_value_for_kernel=i-((kernel_size-1)/2)+j;
	  if (index_value_for_kernel>=0&&index_value_for_kernel<array_size)
	    {
	      if(array[index_value_for_kernel]!=invalid)
		{
		  sum_weight=sum_weight+kernel[j]; // to know the total weigth from kernel that was used
		  sum_value=sum_value+(array[index_value_for_kernel]*kernel[j]); // sum of weigthed value
		}
	    }
	}
      results[i]=sum_value/sum_weight;
    }
  // copy the results to array
  for(int i = 0; i < array_size; i++)
    {
      // if the value was invalid, then leave as it is
      if (array[i]!=invalid)
	{
	  array[i]=results[i];
	}
    }
  free(kernel);
  free(results);
}
SEXP smooth_double_gaussian_cwrap(SEXP array_r, SEXP array_size_r, SEXP sd_r, SEXP invalid_r)
{

  // create the list of variable of correct c types
  double* array;
  int array_size;
  double sd; 
  double invalid;
  /* // coersion */
  array_size=INTEGER_VALUE(array_size_r);
  sd=REAL(sd_r)[0];
  invalid=REAL(invalid_r)[0];

  // by convention, R function do not change the argument, so we will do the same
  array= (double*)malloc(array_size*sizeof(double));
  for(int i=0;i< array_size;i++)
    array[i]=REAL(array_r)[i];

  // call the c function
  smooth_double_gaussian(array, array_size, sd, invalid);

  SEXP out = PROTECT(allocVector(REALSXP, array_size));
  
  for(int i=0;i< array_size;i++)
    REAL(out)[i]=array[i];
  
  UNPROTECT(1);
  return(out);
}

void smooth_double_gaussian_degrees(double* array, int array_size,double smooth, double invalid)
{
  /* Smooth a double array of size "arraysize" using a Gaussian kernal
     Arguments:
     double *array         pointer to data to be smoothed (memory must be pre-allocated)
     int size              number of elements in original array
     int smooth            standard deviation of the kernel
     int invalid           invalid data value to be excluded from smoothing
     
     the gaussian kernel has a size of 3 standard deviation on both side of the middle

     the degrees are changed to radians, then we get the cos and sin of the radians, then smooth the cos and sin, get back the radian with acos and asin and then the degree
  */
  double* s;
  double* c;
  if(array_size<=0)
    return;
  if(smooth<=0)
    return;
  s= (double*)malloc(array_size*sizeof(double));
  c = (double*)malloc(array_size*sizeof(double));
  // check on the input data
  for(int i =0; i < array_size;i++)
    {
      if(array[i]!=-1.0&&array[i]<0)
	{
	  printf(": function smooth_double_gaussian_degrees, there are values that are not -1 and negative\n");
	  return;
	}
      if(array[i]>360)
	{
	  printf("function smooth_double_gaussian_degrees, there are values that are larger than 360\n");
	  return;
	}
    }
    for(int i =0; i < array_size;i++)
    {
      if(array[i]==-1.0)
	{
	  c[i]=-1.0;
	  s[i]=-1.0;
	}
      else
	{
	  c[i]=cos(array[i]/180*M_PI);
	  s[i]=sin(array[i]/180*M_PI);
	}
    }
    // smooth the two linear arrays
    smooth_double_gaussian(c,array_size,smooth,invalid);
    smooth_double_gaussian(s,array_size,smooth,invalid);    
    // get back the degree
    for(int i =0; i < array_size;i++)
    {
      if(array[i]!=-1)
	{
	  if(s[i]>=0)
	    array[i]=acos(c[i])*180/M_PI;
	  if(s[i]<0)
	    array[i]=360-(acos(c[i])*180/M_PI);
	}
    }
    free(c);
    free(s);
    return;
}
SEXP smooth_double_gaussian_degrees_cwrap(SEXP array_r, SEXP array_size_r, SEXP sd_r, SEXP invalid_r)
{
  // create the list of variable of correct c types
  double* array;
  int array_size;
  double sd; 
  double invalid;
  /* // coersion */
  array_size=INTEGER_VALUE(array_size_r);
  sd=REAL(sd_r)[0];
  invalid=REAL(invalid_r)[0];

  // by convention, R function do not change the argument, so we will do the same
  array= (double*)malloc(array_size*sizeof(double));
  for(int i=0;i< array_size;i++)
    array[i]=REAL(array_r)[i];

  // call the c function
  smooth_double_gaussian_degrees(array, array_size, sd, invalid);

  SEXP out = PROTECT(allocVector(REALSXP, array_size));
  
  for(int i=0;i< array_size;i++)
    REAL(out)[i]=array[i];
  
  UNPROTECT(1);
  return(out);
}

SEXP smooth_double_gaussian_2d_cwrap(SEXP array_r, SEXP x_size_r,SEXP y_size_r,SEXP smooth_r,SEXP invalid_r)
{
  double* array;
  double* o;
  int x_size;
  int y_size;
  SEXP out;
  
  array=REAL(array_r);
  x_size=INTEGER_VALUE(x_size_r);
  y_size=INTEGER_VALUE(y_size_r);
  int total=x_size*y_size;
  o=(double*)malloc(total*sizeof(double));

  for(int i = 0;i<total;i++)
    o[i]=array[i];

  smooth_double_gaussian_2d(o,x_size,y_size,REAL(smooth_r)[0],REAL(invalid_r)[0]);
  out = PROTECT(allocMatrix(REALSXP,y_size,x_size));
  double* ans = REAL(out);
  for(int i = 0; i < total; i++)
    ans[i]=o[i];

  free(o);
  UNPROTECT(1);
  return(out);
}


void smooth_double_gaussian_2d(double* array, int x_size,int y_size, double smooth, double invalid)
{
  /* Smooth a 2d array of size "x_size*y_size" using a Gaussian kernal
     Arguments:
        double *array         pointer to data to be smoothed (memory must be pre-allocated)
        int x_size              number of x bins in original array
        int y_size              number of y bins in original array
        int smooth            standard deviation of the kernel
        int invalid           invalid data value to be excluded from smoothing

	the gaussian kernel has a size of 3 standard deviation on both side of the middle
  */
  if(x_size<=0 || y_size<=0)
    {
      printf("in smooth_double_gaussian_2d x_size or y_size of data <=0\n");
      return;
    }
  if(smooth==0)
    {
      return; // do nothing 
    }
  if(smooth<=0)
    {
      printf("in smooth_double_gaussian_2d standard deviation of gaussian kernel is <= 0\n");
      return;
    }

  int num_standard_deviations_in_kernel=3;

  int kernel_size_x=((int)(smooth*num_standard_deviations_in_kernel*2))+1; // for both side
  int kernel_size_y=((int)(smooth*num_standard_deviations_in_kernel*2))+1; // for both side
  if (kernel_size_x%2!=1) // should be an odd number
    { kernel_size_x++; }
  if (kernel_size_y%2!=1) // should be an odd number
    { kernel_size_y++; }
  int kernel_size=kernel_size_x*kernel_size_y;
  double* kernel; // for the gaussian kernel
  double* results; // to store the temporary results 
  double sum_weight;
  double sum_value;
  int x_index_value_for_kernel;
  int y_index_value_for_kernel;
  // make a gaussian kernel
  kernel= (double*)malloc(kernel_size*sizeof(double));
  results =  (double*)malloc(x_size*y_size*sizeof(double));

  gaussian_kernel_2d(kernel,
		     kernel_size_x,
		     kernel_size_y,
		     smooth);// standard deviation in kernel

  // for each x bin
  for(int x = 0; x < x_size; x++)
    {
      // for each y bin
      for (int y = 0; y < y_size; y++)
	{
	  sum_weight=0;
	  sum_value=0;
	  // loop for all the bins in the kernel
	  for (int xx = 0; xx < kernel_size_x; xx++)
	    {
	      for (int yy = 0; yy < kernel_size_y; yy++)
		{                       
		  // find the bin in the data that correspond for that bin in the kernel
		  x_index_value_for_kernel=x-((kernel_size_x-1)/2)+xx;
		  y_index_value_for_kernel=y-((kernel_size_y-1)/2)+yy;
		  // check if that bin is within the original map
		  if (x_index_value_for_kernel>=0 && x_index_value_for_kernel<x_size &&
		      y_index_value_for_kernel>=0 && y_index_value_for_kernel<y_size)
		    {
		      if(array[x_index_value_for_kernel*y_size+y_index_value_for_kernel]!=invalid)
			{
			  sum_weight=sum_weight+kernel[(xx*kernel_size_y)+yy]; //total weigth from kernel that was used
			  sum_value=sum_value+(array[x_index_value_for_kernel*y_size+y_index_value_for_kernel] *
					       kernel[xx*kernel_size_y+yy]); // sum of weigthed value
			}
		    }
		}
	    }
	  // we have looped for the entire kernel, now save the value in results array
	  results[x*y_size+y]=sum_value/sum_weight;
	}
    }
  // copy the results to array
  for(int i = 0; i < x_size*y_size; i++)
    {
      // if the value was invalid, then leave as it is
      if (array[i]!=invalid)
	{
	  array[i]=results[i];
	}
    }
  free(kernel);
  free(results);
}


void gaussian_kernel_2d(double* kernel,
			int x_size,
			int y_size,
			double standard_deviation)
{
  /*function to make a 2d gaussian kernel*/
  if (x_size%2==0)
    {
      printf("in gaussian_kernel_2d, x_size should be an odd number to get a 0 bin\n");
      return;
    }
  if (y_size%2==0)
    {
      printf("in gaussian_kernel, y_size should be an odd number to get a 0 bin\n");
      return;
    }
  if (standard_deviation<=0)
    {
      printf("standard deviation for gaussian_kernel is <= 0\n");
      return;
    }
  int x,y;
  double part_1;
  double part_2;
  double num_part_2;
  double den_part_2;
  part_1=1.0/(2*M_PI*pow(standard_deviation,2)); // should there be 2 pow?
  for (int i = 0; i < x_size; i++)
    {
      for (int j = 0; j < y_size; j++)
	{
	  x=i-(x_size-1)/2; // distance from middle point x
	  y=j-(y_size-1)/2; // distance from middle point y
	  num_part_2=pow(x,2)+pow(y,2);
	  den_part_2=2*pow(standard_deviation,2);
	  part_2= exp(-(num_part_2/den_part_2));
	  kernel[i*y_size+j]=part_1*part_2;
	}
    }
}





void set_array_to_value_int (int* array, int array_size, int value)
{
  /* set all the data of an array to a specific value */
  for (int i = 0; i < array_size; i++)
    {
      array[i]=value;
    }
}
void set_array_to_value_double (double* array, int array_size, double value)
{
  /* set all the data of an array to a specific value */
  for (int i = 0; i < array_size; i++)
    {
      array[i]=value;
    }
}
double correlation (double* x, double* y, int size, double invalid)
{
  /* return the r value of a linear correlation
     see NI Fisher page 145, 6.19 */
  double sum_x=0;
  double sum_y=0;
  double mean_x=0;
  double mean_y=0;
  double sum_x_mean=0;
  double sum_y_mean=0;
  double sum_prod_diff_mean=0;
  double r;
  int n=0;
  for (int i = 0; i < size; i++)
    {
      if (x[i]!=invalid && y[i]!=invalid)
	{
	  sum_x=sum_x+x[i];
	  sum_y=sum_y+y[i];
	  n++;
	}
    }
  mean_x=sum_x/n;
  mean_y=sum_y/n;
  for (int i = 0; i < size; i++)
    {
      if (x[i]!=invalid && y[i]!=invalid)
	{
	  sum_x_mean=sum_x_mean+pow((x[i]-mean_x),2);
	  sum_y_mean=sum_y_mean+pow((y[i]-mean_y),2);
	  sum_prod_diff_mean=sum_prod_diff_mean+((x[i]-mean_x)*(y[i]-mean_y));
	}
    }
  if (sum_x_mean == 0 || sum_y_mean == 0)
    {
      r = 0;
    }
  else
    {
      r=sum_prod_diff_mean/sqrt((sum_x_mean*sum_y_mean));
    }
  if (r<-1||r>1)
    {
      printf("problem with correlation function, value of r out of range\n");
    }
  return r;
}
