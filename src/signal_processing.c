#include <fftw3.h>
#include <math.h>
#include "relectro.h"


/****************************************************************************
 band_pass_filter_one_channel_fftw, for DOUBLE
Take the raw data from one channel, do the fft and then band pass filter 
uses the fftw code
*****************************************************************************/
int band_pass_filter_one_channel_fftw(double* channel_data, int num_samples, int filtered_signal_size,
                                      int sampling_rate,double lower_pass, double higher_pass)
{
  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Because FFTW3 works faster with n is a pow of 2              //
  // we need to do the entire channel in a few steps                //
  // use windows that overlap                                                   //
  // **********                                                   //
  //        **********                                            //
  //               ***********                                    //
  //                       ***********                            //
  //                                                              //
  // the overlap prevent the artefacts that occurs at the end of  //
  // the filtered segment                                         //
  //                                                              //
  // or if the segment is too short, then pad with 0              //
  //////////////////////////////////////////////////////////////////
  int i, j; // counters
  // checks for arguments
  if (num_samples <=0)
  {
    fprintf(stderr,"band_pass_filter_one_channel_fftw\n");
    fprintf(stderr,"num_samples should be larger than 0\n");
    return -1;
  }
  if (sampling_rate <=0 || sampling_rate>=100000)
  {
    Rprintf("band_pass_filter_one_channel_fftw\n");
    Rprintf("sampling_rate should be between 0 and 100000 Hz\n You gave %d\n",sampling_rate);
    return -1;
  }
  // should be a power of 2
  if (!(((filtered_signal_size | (filtered_signal_size - 1)) + 1) / 2 == filtered_signal_size))
  {
    Rprintf("band_pass_filter_one_channel_fftw\n");
    Rprintf("filtered_signal_size should be a power of 2\nYou gave %d\n",filtered_signal_size);
    Rprintf("try one of the following values\n");
    for (i = 10; i < 18; i++)
    {
      Rprintf("%d\n",(int)pow((double)2,(double)i));
    }	
    return -1; 
  }
  if(lower_pass <=0 || lower_pass>=sampling_rate/2)
  {
    // if lower_pass is set at 0, the filtered signal will be at 0
    Rprintf("band_pass_filter_one_channel_fftw\n");
    Rprintf("lower_pass needs to be between 0 and %d Hz\n You gave %lf\n",sampling_rate/2,lower_pass);
    return -1;
  }
  if(higher_pass > sampling_rate/2)
  {
    Rprintf("%s: band_pass_filter_one_channel_fftw\n");
    Rprintf("higher_pass is larger than sampling_rate/2\n You gave %lf\n",higher_pass);
    return -1;
  }
  if(lower_pass >= higher_pass)
  {
    Rprintf("band_pass_filter_one_channel_fftw\n");
    Rprintf("lower_pass should be smaller than higher_pass\n You gave %lf as min and %lf as max\n",lower_pass,higher_pass);
    return -1;
  }
  // in fftw3, m has a different meaning than NR3
  // in fftw3 it is the size of the complex array returned by 
  // the fft
  int m = filtered_signal_size/2+1;
  int overlap; // size of overlap between two windows at one end
  if (filtered_signal_size < 10000)
  {
    overlap=filtered_signal_size/2; // to minimise artefacts
  }
  else
  {
    overlap=filtered_signal_size/3; 
  }
  int non_overlap=filtered_signal_size-overlap; // the beginning of window until overlap
  int window_end_chopped=overlap/2; //what will not be in the data in the end 
  int num_windows;
  double* filter_function; // kernel to do the filtering after the fft
  filter_function = (double*) malloc(sizeof(double)*m);
  fftw_complex *out; // complex array return by fft forward
  double* filtered_signal; // will go in and out of the fft
  fftw_plan fft_plan_forward; // plan to do fft forward,
  fftw_plan fft_plan_backward; 
  out= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * filtered_signal_size);
  filtered_signal=(double*) fftw_malloc(sizeof(double) * filtered_signal_size);
  // create the plan to do the fft
  fft_plan_forward= fftw_plan_dft_r2c_1d(filtered_signal_size, filtered_signal, out,FFTW_MEASURE);
  fft_plan_backward= fftw_plan_dft_c2r_1d(filtered_signal_size, out, filtered_signal,FFTW_MEASURE);
  make_butterworth_filter(sampling_rate, // max freq is sr/2
                          m, // m
                          filter_function, // size m
                          lower_pass, // lower cut-off frequency
                          higher_pass); // higher cut-off frequency
  int start, end; // index for start and end of window
  if (num_samples>filtered_signal_size)// there is more than one window for the fft
  {
    // calculate the number of windows needed to fft channel //
    end=filtered_signal_size-window_end_chopped; // first window has that size
    num_windows=1;
    while (end < num_samples)
    {
      end=end+non_overlap; // the subsequent windows are smaller because cut at both ends
      num_windows++;
    }
    // loop and do fft and filtering for each window and then copy in 
    // channel data
    for (i = 0; i < num_windows; i++)
    {
      if (i==0)
      {
        start=0;
        end=filtered_signal_size-window_end_chopped;
      }
      if (i != num_windows-1&&i!=0)
      {// not the last, not the first
        start = (i * non_overlap)-window_end_chopped; // first window is window_end_chopped smaller than other
        end = start + filtered_signal_size;
      }
      if (i== num_windows-1) // last block
      {
        // needs to be correct size but ends at the end of data
        end = num_samples - 1;
        start = end - filtered_signal_size;
      }
      // put the data into the vector
      for (j = 0; j < filtered_signal_size; j++)
      {
        filtered_signal[j] =  channel_data[start+j];
      }
      // do the fft forward, will give us an array of complex values
      fftw_execute(fft_plan_forward); // results is stored in out array
      
      // do the filtering, array out has a size of m
      for(j=0; j < m;j++)
      {
        out[j][0]=out[j][0]*filter_function[j];
        out[j][1]=out[j][1]*filter_function[j];
      }
      // reverse the fft
      fftw_execute(fft_plan_backward);
      // rescale the results by the factor of filtered_signal_size
      for (j=0;j < filtered_signal_size;j++)
      {
        filtered_signal[j]=filtered_signal[j]/filtered_signal_size;
      }
      // copy the filtered data back to channel_data array
      if (i==0 && i != num_windows-1) // first window but not last window
      {
        for (j = 0; j < filtered_signal_size-window_end_chopped; j++)
        {
          // from beginning but not the end that will be in the next segment
          channel_data[start+j] = (short)filtered_signal[j];
        }
      }
      if (i!=0 && i != num_windows-1) // middle window
      {
        for (j = window_end_chopped; j < filtered_signal_size - window_end_chopped; j++)
        {
          channel_data[start+j] = (short)filtered_signal[j];
        }
      }
      if (i==num_windows-1) // last window
      {
        for (j = window_end_chopped; j <= filtered_signal_size; j++)
        {
          channel_data[start+j] = (short)filtered_signal[j];
        }
      } 
    }
  }
  /// that means that we just do one shot padded with 0
  if (filtered_signal_size>num_samples)
  {
    // read the data and then padd with 00000
    for (j = 0; j < num_samples; j++)
    {
      filtered_signal[j] = channel_data[j];
    }
    for (j = num_samples; j < filtered_signal_size; j++)
    {
      filtered_signal[j] =  0;
    }
    // do the fft forward, will give us an array of complex values
    fftw_execute(fft_plan_forward); // results is stored in out array
    // do the filtering, array out has a size of m
    for(j=0; j < m;j++)
    {
      out[j][0]=out[j][0]*filter_function[j];
      out[j][1]=out[j][1]*filter_function[j];
    }
    // reverse the fft
    fftw_execute(fft_plan_backward);
    // rescale the results by the factor of filtered_signal_size
    for (j=0;j < filtered_signal_size;j++)
    {
      filtered_signal[j]=filtered_signal[j]/filtered_signal_size;
    }
    for (j=0; j < num_samples; j++)
    {
      channel_data[j]=(short)filtered_signal[j];
    }
  }
  free(filter_function);
  fftw_destroy_plan(fft_plan_forward);
  fftw_destroy_plan(fft_plan_backward);
  fftw_free(filtered_signal); fftw_free(out);
  return 0;
}


/*****************************************************************
/// function to make a filter function for a band pass filter ////
/// a low and high pass filter function is made and then      ////
/// they are multiplied together to give the band pass filter ////
/// function                                                  ////
*****************************************************************/
int make_butterworth_filter(int sampling_rate, // max freq is sr/2
                             int filter_length, //
                             double* filter_function, // size m
                             double low_pass, // lower cut-off frequency
                             double high_pass) // higher cut-off frequency
{
  int n_low; // order of filter
  int n_high; // order of filter
  double* function_low_pass; 
  double* function_high_pass;
  double frequency_steps=(double)sampling_rate/2/(double)filter_length;
  double frequency;
  int size_double=sizeof(double);
  if (low_pass < 0 || low_pass > sampling_rate/2)
  {
    Rprintf("lower_pass should be between 0 and %d\n", sampling_rate/2 );
    Rprintf("make_butterworth_filter()\n");
    return 1;
  }
  if (high_pass < low_pass || high_pass > sampling_rate/2)
  {
    Rprintf("high_pass should be between %d and %d\n", low_pass, sampling_rate/2);
    Rprintf("is %d\n",high_pass);
    Rprintf("make_butterworth_filter()\n");
    return 1;
  }
  // allocate memory
  function_low_pass = (double*)malloc(size_double*filter_length); // returns pointer to first element
  function_high_pass = (double*)malloc(size_double*filter_length);
  // set the order of the filters
  if (low_pass>=0&&low_pass<15)
  {
    n_low=15;
  }
  if (low_pass>=15 &&low_pass<85)
  {
    n_low=19;
  }
  if (low_pass>=85)
  {
    n_low=19;
  }
  if (high_pass>=0&&high_pass<15)
  {
    n_high=15;
  }
  if (high_pass>=15 &&high_pass<85)
  {
    n_high=19;
  }
  if (high_pass>=85)
  {
    n_high=19;
  }
  // make the function for low pass
  frequency=0;
  for (int i = 0; i < filter_length; i++)
  {      
    function_low_pass[i]= 1/sqrt(1+pow((frequency/high_pass),2*n_high));
    frequency=frequency+frequency_steps;
  }
  // make the function for high pass
  frequency=0;
  for (int i = 0; i < filter_length; i++)
  {
    function_high_pass[i]=1-(1/sqrt(1+pow((frequency/low_pass),2*n_low)));
    frequency=frequency+frequency_steps;
  }
  // make the product of the two function
  frequency=0;
  for (int i = 0; i < filter_length; i++)
  {
    filter_function[i]=function_high_pass[i]*function_low_pass[i];
  }
  free(function_low_pass);
  free(function_high_pass);
  return 0;
}

SEXP band_pass_filter_one_channel_fftw_cwrap(SEXP channel_data_r, SEXP num_samples_r,
                                            SEXP sampling_rate_r, SEXP lower_pass_r, SEXP higher_pass_r){
  int num_samples=INTEGER_VALUE(num_samples_r);
  int sampling_rate = INTEGER_VALUE(sampling_rate_r);
  double lower_pass = REAL(lower_pass_r)[0];
  double higher_pass = REAL(higher_pass_r)[0];
  int filtered_signal_size=2;
  
  if(num_samples<100000){
    while(filtered_signal_size<num_samples)
      filtered_signal_size=filtered_signal_size*2;
  } 
  else {
    filtered_signal_size=pow(2,17);
  }
  
  // will not affect the input vector
  double* array= (double*)malloc(num_samples*sizeof(double));
  for(int i=0;i< num_samples;i++)
    array[i]=REAL(channel_data_r)[i];
  
  band_pass_filter_one_channel_fftw(array, num_samples, 
                                    filtered_signal_size,
                                    sampling_rate, lower_pass, higher_pass);
  
  SEXP out = PROTECT(allocVector(REALSXP, num_samples));
  for(int i=0;i< num_samples;i++)
    REAL(out)[i]=array[i];
  UNPROTECT(1);
  free(array);
  return out;
}

SEXP power_root_mean_square(SEXP channel_data_r, SEXP num_samples_r, 
                            SEXP window_size_samples_r, SEXP window_slide_r){
  int num_samples=INTEGER_VALUE(num_samples_r);
  int window_size_samples=INTEGER_VALUE(window_size_samples_r);
  int window_slide=INTEGER_VALUE(window_slide_r);
  
  if(num_samples<0){
    Rprintf("power_root_mean_square: num_samples should be larger than 0%d\n", num_samples);
    return 0;
    
  }
  if(window_size_samples<0){
     Rprintf("power_root_mean_square: window_size_samples should be larger than 0%d\n", window_size_samples);
    return 0;
  }
  if(window_slide<0){
    Rprintf("power_root_mean_square: window_slide should be larger than 0%d\n", window_slide);
    return 0;
  }
  
  int num_windows=1+(num_samples-window_size_samples)/window_slide;
  
  double* pwr= (double*)malloc(num_windows*sizeof(double));
  int start_index;
  int end_index;
  double sum_square;
  double mean_square;
  
  for(int i = 0; i < num_windows; i++){
    start_index=i*window_slide;
    end_index=start_index+window_size_samples;
    sum_square = 0;
    
    for (int j = start_index; j < end_index; j++)
    {
      sum_square=sum_square+(REAL(channel_data_r)[j]*REAL(channel_data_r)[j]);
    }
    mean_square=sum_square/window_size_samples;
    pwr[i]=sqrt(mean_square);
  }
  
  SEXP out = PROTECT(allocVector(REALSXP, num_windows));
  for(int i=0;i< num_windows;i++)
    REAL(out)[i]=pwr[i];
  free(pwr);
  UNPROTECT(1);
  return(out);
}
