#include <fftw3.h>
#include <math.h>
#include "relectro.h"

void convolution_fftw3(double* signal, int signal_size, double* response, int response_size, double* out)
{
  int output_size=signal_size; // size of the vector out
  if(output_size<response_size){
    printf("convolution_fftw3, output_size set to response_size\n");
    output_size=response_size;
  }
  
  //printf("convolution_fftw3\n");
  //printf("signal_size: %d\n",signal_size);
  //printf("response_size: %d\n",response_size);
  //printf("output_size: %d\n",output_size);

  double* signal_in;
  double* response_in;
  
  // to ensure there is not wrapping in the convolution
  int in_size1 = (signal_size+response_size)-1;
  int in_size2 = in_size1;
  // get a power of 2 
  int k =1;
  while(in_size2 > pow((double)2,(double)k))
    k++;
  in_size2=pow((double)2,(double)k);
 // printf("convolution in_size1: %d\n",in_size1);
//  printf("convolution in_size2: %d\n",in_size2);
  signal_in = (double*) fftw_malloc(sizeof(double) * in_size2);
  response_in= (double*) fftw_malloc(sizeof(double) * in_size2);
  
  // get the signal with padding
  for(int i = 0; i < signal_size;i++)
    signal_in[i]=signal[i];
  for(int i = signal_size; i < in_size2;i++)
    signal_in[i]=0;
  // get the response with padding
  for(int i = 0; i < response_size;i++)
    response_in[i]=response[i];
  for(int i = response_size; i < in_size2;i++)
    response_in[i]=0;
  
  fftw_complex *signal_out;
  fftw_complex *response_out;
  fftw_complex *convolution_in;
  
  signal_out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (in_size2/2+1)); // output of r2c n/2+1
  response_out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (in_size2/2+1)); // output of r2c n/2+1
  convolution_in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (in_size2/2+1)); // output of r2c n/2+1
  
  fftw_plan plan_signal;
  plan_signal=fftw_plan_dft_r2c_1d(in_size2, signal_in, signal_out, FFTW_ESTIMATE);
  fftw_plan plan_response;
  plan_response=fftw_plan_dft_r2c_1d(in_size2, response_in, response_out, FFTW_ESTIMATE);
  
  fftw_execute(plan_signal);
  fftw_execute(plan_response);

  // multiply the complex signal and response
  int m = in_size2/2+1;
  double scale=1.0/(double)in_size2;
  for(int i = 0; i < m;i++)
  {
    convolution_in[i][0] = (signal_out[i][0] * response_out[i][0]
                           - signal_out[i][1] * response_out[i][1])* scale;
    convolution_in[i][1] = (signal_out[i][0] * response_out[i][1]
                           + signal_out[i][1] * response_out[i][0])* scale;
  }
  
  fftw_plan plan_convolution;
  plan_convolution = fftw_plan_dft_c2r_1d(in_size2, convolution_in, signal_in, FFTW_ESTIMATE);
  
  // fftw inverse the results
  fftw_execute(plan_convolution);

  // remove the padding to get to the next pow2
  int power_pad=in_size2-in_size1;
  for(int i = 0; i < output_size ;i++)
  {
    out[i]=signal_in[i+response_size/2]; // Not sure if response_size/2 is ok for offset
  }
  
  fftw_destroy_plan(plan_signal);
  fftw_destroy_plan(plan_response);
  fftw_destroy_plan(plan_convolution);
  fftw_free(signal_in);
  fftw_free(response_in);
  fftw_free(signal_out);
  fftw_free(response_out);
  fftw_free(convolution_in);
}

SEXP convolution_fftw3_cwrap(SEXP inA_r,SEXP inA_size_r, SEXP inB_r, SEXP inB_size_r)
{
  int output_size=INTEGER_VALUE(inA_size_r);
  if(output_size<INTEGER_VALUE(inB_size_r))
    output_size=INTEGER_VALUE(inB_size_r);

  printf("convolution_fftw3_cwrap\n");  
  SEXP out = PROTECT(allocVector(REALSXP, output_size));
  double* oo = REAL(out);
  convolution_fftw3(REAL(inA_r),INTEGER_VALUE(inA_size_r),REAL(inB_r),INTEGER_VALUE(inB_size_r),oo);
    
  UNPROTECT(1);
  return(out);
}

			    
