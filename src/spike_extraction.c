#include <math.h>
#include "relectro.h"

SEXP identify_spike_times(SEXP dataf_r,
                          SEXP dataf_size_r,
                          SEXP power_r,
                          SEXP power_size_r,
                          SEXP powerWindowSize_r,
                          SEXP powerWindowSlide_r,
                          SEXP powerThreshold_r)
  
  {
  
  int dataf_size=INTEGER_VALUE(dataf_size_r);
  int power_size=INTEGER_VALUE(power_size_r);
  int powerWindowSize=INTEGER_VALUE(powerWindowSize_r);
  int powerWindowSlide=INTEGER_VALUE(powerWindowSlide_r);
  double powerThreshold=REAL(powerThreshold_r)[0];
  double* power = REAL(power_r);
  double* dataf = REAL(dataf_r);
  
  if(dataf_size<0){
    Rprintf("identify_spike_time: dataf_size <0, %d\n", dataf_size);
    return 0;
  }
  if(power_size<0){
    Rprintf("identify_spike_time: power_size <0, %d\n", power_size);
    return 0;
  }
  
  int num_windows=1+(dataf_size-powerWindowSize)/powerWindowSlide;
  if(num_windows!=power_size){
    Rprintf("identify_spike_time: power_size = %d but should be %d\n", power_size,num_windows);
    return 0;
  }

  
  
  // set the arry to one spike per power window, which is way too many
  double* negVal = (double*)malloc(num_windows*sizeof(double));
  int* spikes=(int*)malloc(num_windows*sizeof(int));
  double* powerPeak = (double*)malloc(num_windows*sizeof(double));
  int numSpikes=0;
  
  
  // loop in the power array
  int j;
  int indexDataStart;
  int indexDataEnd;
  for(int i = 0; i < num_windows; i++)
  {
    if(power[i]>powerThreshold){
      j=i;
      while(i<num_windows&&power[i]>powerThreshold) //loop until power is below threshold
        i++;
      i--; // don't consider the last while iteration that was below threshold
      // above threshold from index j to i
      
      // get the peak power
      for(int k = j; k <= i; k++)
      {
        if(k==j)
          powerPeak[numSpikes]=power[k];
        else{
          if(powerPeak[numSpikes]<power[k])
            powerPeak[numSpikes]=power[k];
        }
      }
      
      
      indexDataStart=powerWindowSlide*j;
      indexDataEnd=powerWindowSlide*i+powerWindowSize;
      for(int k=indexDataStart; k < indexDataEnd; k++)
      {
        if(k==indexDataStart){
          spikes[numSpikes]=k;
          negVal[numSpikes]=dataf[k];
        }
        else{
          if(dataf[k]<negVal[numSpikes])
          {
            spikes[numSpikes]=k;
            negVal[numSpikes]=dataf[k];
          }
        }
      }
      numSpikes++;
    }
  }
  
  SEXP out = PROTECT(allocMatrix(REALSXP,numSpikes,3));
  double* outc = REAL(out);
  for(int i = 0; i < numSpikes; i++){
    outc[i]=spikes[i];
    outc[i+numSpikes]=negVal[i];
    outc[i+numSpikes*2]=powerPeak[i];
  }
  free(spikes);
  free(negVal);
  free(powerPeak);
  UNPROTECT(1);
  return(out);
}
