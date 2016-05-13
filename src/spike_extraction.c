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
      
      // get the peak power, not necessarily at the same time as spike
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





SEXP merge_simultaneous_spikes(SEXP time_r,
                               SEXP trough_r,
                               SEXP size_r,
                               SEXP max_time_difference_r)
{
  // merge spikes with isi smaller than max_time_difference
  // when merging use spike time of the spike with smallest trough
  int size=INTEGER_VALUE(size_r);
  int max_time_difference=INTEGER_VALUE(max_time_difference_r);
  int* time = INTEGER_POINTER(time_r);
  double* trough = REAL(trough_r);
  
  if(size<=0){
    Rprintf("merge_simultaneous_spike: size <=0, %d\n", size);
    return 0;
  }
  
  // set the arry to one spike per power window, which is way too many
  int* ti=(int*)malloc(size*sizeof(int));
  int num_spikes=0;
  
  int j;
  int index_keep=0;
  int smallest_trough=0;
  for(int i = 0; i < size; i++)
  {
   // printf("i:%d, time[i]: %d\n",i,time[i]);
    j=i;
    while((j+1)<size && (time[j+1]-time[i])<max_time_difference){
      //printf("join with j:%d, time[j]: %d\n",j,time[j]);
      j++;
    }
    
  //  printf("from %d to %d\n",i,j);
    
    if(i==j){ // no simultaneous spikes
      ti[num_spikes]=time[i];
      //printf("adding single i == %d,  num_spikes: %d, ti[num_spikes]: %d\n",i,num_spikes,ti[num_spikes]);
      num_spikes++;
    }
    else{ // simultaneous spikes (from index i to j), find the spike with smallest trough and keep it.
      for(int k = i; k <= j; k++){
        if(k==i){
          index_keep=k;
          smallest_trough=trough[k];
        }
        else{
          if(trough[k]<smallest_trough){
            index_keep=k;
            smallest_trough=trough[k];
            }
          }
        }
        ti[num_spikes]=time[index_keep];
      //  printf("adding single for i:%d to j:%d, num_spikes: %d, ti[num_spikes]: %d\n",i,j,num_spikes,ti[num_spikes]);
        num_spikes++;
      }
    i=j;
    }
    
  SEXP out = PROTECT(allocVector(REALSXP, num_spikes));
  double* outc = REAL(out);
  for(int i = 0; i < num_spikes; i++)
    outc[i]=ti[i];

  free(ti);
  UNPROTECT(1);
  return(out);
}

SEXP create_spk_file(SEXP data_r, SEXP nrow_r, SEXP ncol_r, SEXP res_r, SEXP res_lines_r, SEXP window_r, SEXP file_name_r){
  const char* file_name = CHAR(STRING_ELT(file_name_r,0));
  int * m = INTEGER_POINTER(data_r);
  int nrow=INTEGER_VALUE(nrow_r);
  int ncol=INTEGER_VALUE(ncol_r);
  int* res=INTEGER_POINTER(res_r);
  int res_lines = INTEGER_VALUE(res_lines_r);
  int window = INTEGER_VALUE(window_r);
  
  
  Rprintf("writing %s from %d %d array and %d spikes\n",
          file_name,
          nrow,ncol,
          res_lines);
  
  // memory for one spike
  short int* out=(short int*)malloc(sizeof(short int)*window*ncol);
  
  FILE *my_file = fopen(file_name, "wb");
  
  //
  int index;
  int wb;// begin window
  for(int i = 0; i < res_lines; i++){
    index=0;
    wb=res[i]-window/2;
    for(int k =0; k < window; k++){
      for(int j = 0; j < ncol; j++){
        out[index]= (short int)m[j*nrow+wb+k];
        index++;
      }
    }
   fwrite(out,sizeof(short int),window*ncol,my_file);
  // save this spike
  }
  fclose(my_file);
  free(out);  
  return(R_NilValue);
}


SEXP get_waveform_matrix(SEXP signal_r, SEXP signal_lines_r, SEXP res_r, SEXP res_lines_r, SEXP window_r){
  int * signal = INTEGER_POINTER(signal_r);
  int signal_lines =INTEGER_VALUE(signal_lines_r);
  int* res=INTEGER_POINTER(res_r);
  int res_lines = INTEGER_VALUE(res_lines_r);
  int window = INTEGER_VALUE(window_r);
  
  SEXP out = PROTECT(allocMatrix(INTSXP,res_lines,window));
  int* o = INTEGER_POINTER(out);
  
  int wb;// begin window
  for(int i = 0; i < res_lines; i++){
    wb=res[i]-window/2;
    for(int k =0; k < window; k++){
        o[k*res_lines+i]= signal[wb+k];
      }
    }

  UNPROTECT(1);
  return(out);
}