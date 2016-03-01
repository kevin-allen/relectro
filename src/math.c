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
