#include <math.h>
#include "relectro.h"

SEXP identify_spike_times(SEXP dataf_r,
                          SEXP dataf_size_r,
                          SEXP power_r,
                          SEXP power_size_r,
                          SEXP powerWindowSize_r,
                          SEXP powerWindowSlide_r,
                          SEXP powerThreshold_r,
                          SEXP refractory_r)
  {
// function to get the spike times from the high power events.
// if the size of the powerWindow is small the trough of the spike
// might fall outside of the powerWindow. This is why, I modify this
// so that the trough can also be just before and after the window with
// high power. Power detect either raising or descending peak of > 700 Hz
// oscillation, so the sholders should be approx 0.5 ms  (~ 10 samples)
// will try with sholders of the size of the refractory period
  

  int dataf_size=INTEGER_VALUE(dataf_size_r);
  int power_size=INTEGER_VALUE(power_size_r);
  int powerWindowSize=INTEGER_VALUE(powerWindowSize_r);
  int powerWindowSlide=INTEGER_VALUE(powerWindowSlide_r);
  int refractory=INTEGER_VALUE(refractory_r);
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
  if(refractory<0){
    Rprintf("identify_spike_time: refractory < 0, %d\n",refractory);
    return(0);
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
  
  int j;
  int indexDataStart;
  int indexDataEnd;
  
  int startWindow; // because of sholder to the power window
  int endWindow; // because of solder to the power window
  j=0;
  while(powerWindowSize*j<=refractory){
    j++;
  }
  startWindow=j;
  endWindow=num_windows-j;
  
  
  
  for(int i = startWindow; i < endWindow; i++)
  {
    if(power[i]>powerThreshold){ // this mean a spike
      j=i;
      while(j+1<num_windows&&power[j+1]>powerThreshold) //loop until power is below threshold
        j++;
      // power is up from i to j, there will only be one spike from i to j
      // get the peak power, not necessarily at the same time as spike
      for(int k = i; k <= j; k++)
      {
        if(k==i)
          powerPeak[numSpikes]=power[k];
        else{
          if(powerPeak[numSpikes]<power[k])
            powerPeak[numSpikes]=power[k];
        }
      }
      
      // get the trace from begining of first high
      indexDataStart=powerWindowSlide*i-refractory;
      indexDataEnd=powerWindowSlide*j+powerWindowSize+refractory; // 0 to 9, 10 to 19
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
      i=j; // will be incremented again with the for loop
      numSpikes++;
    }
  }
  
  
  // make sure that the same spike was not detected twice, and that refactory period was respected
  int* toRemove = (int*)malloc(numSpikes*sizeof(int));
  for(int i =0; i < numSpikes;i++)
    toRemove[i]=0;
  int refUntil;
  for(int i = 0; i < numSpikes-1;i++) // there should be at least one spike after current spike
  {
    if(toRemove[i]==0){ // we have not already decided to remove this spike
      refUntil=spikes[i]+refractory; // do not accept spike until refUntil
      j=i+1;// move to next spike
      while(j<numSpikes&&spikes[j]<refUntil)
      {
        //Rprintf("%d %d %d %d %d\n",i,j,spikes[i], spikes[j], refUntil);
        toRemove[j]=1;
        j++;
      }
    }
  }
  int numValidSpikes=0;
  for(int i = 0; i < numSpikes;i++)
  {
    if(toRemove[i]==0){
      spikes[numValidSpikes]=spikes[i];
      negVal[numValidSpikes]=negVal[i];
      powerPeak[numValidSpikes]=powerPeak[i];
      numValidSpikes++;
    }
  }
  numSpikes=numValidSpikes;
  
  
  
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
  free(toRemove);
  UNPROTECT(1);
  return(out);
}








SEXP spike_waveform_from_traces(SEXP data_r, SEXP nrow_r, SEXP ncol_r, SEXP res_r, SEXP res_lines_r, SEXP window_r){
  int * m = INTEGER_POINTER(data_r);
  int nrow=INTEGER_VALUE(nrow_r);
  int ncol=INTEGER_VALUE(ncol_r);
  int* res=INTEGER_POINTER(res_r);
  int res_lines = INTEGER_VALUE(res_lines_r);
  int window = INTEGER_VALUE(window_r);
  
  // memory for one spike
  int* o=(int*)malloc(sizeof(int)*window*ncol*res_lines);
  
  //
  int index=0;
  int wb;// begin window
  for(int j = 0; j < ncol; j++){ // channels
    for(int k = 0; k < window; k++){
      for(int i = 0; i < res_lines; i++){ // spikes
          wb=res[i]-window/2;
          o[index]=m[j*nrow+wb+k];
          index++;
      }
    }
  }
  
  int size=window*ncol*res_lines;
  SEXP out = PROTECT(allocVector(INTSXP,size));
  int* outc = INTEGER_POINTER(out);
  for(int i = 0; i < size; i++)
    outc[i]=o[i];
  
  free(o);  
  UNPROTECT(1);
  return(out);
}

SEXP create_spk_file(SEXP data_r, SEXP nrow_r, SEXP ncol_r, SEXP res_r, SEXP res_lines_r, SEXP window_r, SEXP file_name_r,SEXP append_r){
  const char* file_name = CHAR(STRING_ELT(file_name_r,0));
  int * m = INTEGER_POINTER(data_r);
  int nrow=INTEGER_VALUE(nrow_r);
  int ncol=INTEGER_VALUE(ncol_r);
  int* res=INTEGER_POINTER(res_r);
  int res_lines = INTEGER_VALUE(res_lines_r);
  int window = INTEGER_VALUE(window_r);
  int append=INTEGER_VALUE(append_r);
  
  // memory for one spike
  short int* out=(short int*)malloc(sizeof(short int)*window*ncol);
  
  FILE *my_file;
  if(append==1){
   my_file= fopen(file_name, "a");
  }
  else{
    my_file=fopen(file_name, "w");
  }
  
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
  //int signal_lines =INTEGER_VALUE(signal_lines_r);
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


SEXP write_fet_file(SEXP nFeatures_r, SEXP nSpikes_r, SEXP fet_r, SEXP fileName_r, SEXP append_r){
  int nFeatures = INTEGER_VALUE(nFeatures_r);
  int nSpikes = INTEGER_VALUE(nSpikes_r);
  int * fet = INTEGER_POINTER(fet_r);
  const char* file_name = CHAR(STRING_ELT(fileName_r,0));
  int append=INTEGER_VALUE(append_r);
  
  FILE *f;
  if(append==1){
    f = fopen(file_name, "a");
  }else{
    f = fopen(file_name, "w");
  }
  if (f==NULL){
    Rprintf("Error opening %s\n",file_name);
    return(R_NilValue);
  }
  
  fprintf(f, "%d\n",nFeatures);
  for(int i=0; i < nSpikes; i++){
    for(int j=0; j < nFeatures-1; j++){
      fprintf(f,"%d ",fet[j*nSpikes+i]);
    }
    fprintf(f,"%d\n",fet[(nFeatures-1)*nSpikes+i]);
  }
  fclose(f);
  return(R_NilValue);
}

SEXP spike_geometrical_features(SEXP spikes_r,SEXP nSpikes_r, SEXP spikeSize_r){
  double * spikes = REAL(spikes_r);
  int nSpikes = INTEGER_VALUE(nSpikes_r);
  int spikeSize = INTEGER_VALUE(spikeSize_r);
  double * spk; // pointer to single spike
  int featuresPerSpike=5;
  SEXP out = PROTECT(allocMatrix(REALSXP,nSpikes,featuresPerSpike));
  double* o=REAL(out);
  
  double baseline; // mean of the first and last quarter of the waveform
  double trough; // most negative value
  double b1,b2;
  double amplitude;
  double threshold_amplitude;
  double threshold_amplitude_proportion=0.5;
  int trough_index;
  int index;
  double interpolation_start;
  double interpolation_end;
  int start_spike_index;
  int end_spike_index;
  double spike_duration;
  double first_half_duration;
  double second_half_duration;
  double max;
 // int index_peak_from_baseline_pre;
//  int index_peak_from_baseline_post;
  double peak_from_baseline_pre;
  double peak_from_baseline_post;
  double peak_amplitude_asymmetry;
  int upsideDown;
  for(int i = 0 ; i < nSpikes; i++){
    spk = spikes+(i*spikeSize); // point to a single spike
    
    // get baseline voltage
    b1=0; // first baseline
    b2=0; // second baseline
    b1=mean_double(spikeSize/4,spk);
    b2=mean_double(spikeSize/4,spk+spikeSize-(spikeSize/4));
    baseline=(b1+b2)/2;
    
    // get trough
    // we assume that the trough is near the middle of the waveform
    upsideDown=0;
    trough=baseline;
    trough_index=spikeSize/3;
    for(int j = spikeSize/3; j < spikeSize-(spikeSize/3);j++){
      if(trough>spk[j]){
        trough=spk[j];
        trough_index=j;
      }
    }
    if(trough==baseline){ // spike is upside down
      upsideDown=1;
      trough_index=spikeSize/3;
      for(int j = spikeSize/3; j < spikeSize-(spikeSize/3);j++){
        if(trough<spk[j]){
          trough=spk[j];
          trough_index=j;
        }
      } 
    }
    amplitude=baseline-trough;
    
    if(upsideDown==0){
      // get width 
      // need the to know the voltage at which we want to calculate width
      threshold_amplitude=baseline-(amplitude*threshold_amplitude_proportion);
      // find the index of threshold before trough
      index = trough_index;
      interpolation_start = 0; // proportion of a bin to add if we interpolate to the threshold
      interpolation_end = 0;
      while(index>0 && spk[index]<threshold_amplitude)
        index--;
      start_spike_index=index;
      if(index>0)
      {
        interpolation_start= (spk[index]-threshold_amplitude)/(spk[index]-spk[index+1]); // proportion of a bin to add
        interpolation_start=start_spike_index+interpolation_start;
      }else
      {
        interpolation_start=start_spike_index;
      }
      // find the index of threshold after trough
      index = trough_index;
      while(index<spikeSize && spk[index]<threshold_amplitude)
        index++;
      end_spike_index=index;
      if(index<spikeSize-1)
      {
        interpolation_end=(spk[index]-threshold_amplitude)/(spk[index]-spk[index-1]); // proportion of a bin to add
        interpolation_end=end_spike_index-interpolation_end;
      }else{
        interpolation_end=end_spike_index;
      }
    }else{
      // get width 
      // need the to know the voltage at which we want to calculate width
      threshold_amplitude=baseline-(amplitude*threshold_amplitude_proportion);
      // find the index of threshold before trough
      index = trough_index;
      interpolation_start = 0; // proportion of a bin to add if we interpolate to the threshold
      interpolation_end = 0;
      while(index>0 && spk[index]>threshold_amplitude)
        index--;
      start_spike_index=index;
      if(index>0)
      {
        interpolation_start= (spk[index]-threshold_amplitude)/(spk[index]-spk[index+1]); // proportion of a bin to add
        interpolation_start=start_spike_index+interpolation_start;
      }else{
        interpolation_start=start_spike_index;
      }
      // find the index of threshold after trough
      index = trough_index;
      while(index<spikeSize && spk[index]>threshold_amplitude)
        index++;
      end_spike_index=index;
      if(index<spikeSize-1)
      {
        interpolation_end=(spk[index]-threshold_amplitude)/(spk[index]-spk[index-1]); // proportion of a bin to add
        interpolation_end=end_spike_index-interpolation_end;
      }else{
        interpolation_end=end_spike_index;
      }
    }
    
    
    spike_duration=(interpolation_end-interpolation_start);
    first_half_duration=(trough_index-interpolation_start);
    second_half_duration=(interpolation_end-trough_index);
    
    if(spike_duration<=0)
    {
      Rprintf("\n");
      for(int j = 0; j < spikeSize; j++)
        Rprintf("%lf\n",j,spk[j]);
      Rprintf("b1:%lf b2:%lf base:%lf trou:%lf trough_index: %d amp:%lf spikeD:%lf firstHalfD:%lf secondHalfD:%lf\n",
              b1,b2,baseline,trough,trough_index,amplitude,spike_duration,first_half_duration,second_half_duration);
      Rprintf("inter_start:%d inter_end:%d\n",
              interpolation_start,interpolation_end);
      
              
    }
    
    
    
    
    // asymmetry
    // calculate the peak amplitude asymmetry (a-b)/(|a|+|b|)
    // a and b being peak_from_baseline_pre and peak_from_baseline_post
    if(upsideDown==0){
      max=trough;
      for (int i =0; i < trough_index;i++)
      {
        if(spk[i]>max)
        {
          max=spk[i];
          // index_peak_from_baseline_pre=i;
        }
      }
      peak_from_baseline_pre=max-baseline;
      max=trough;
      for (int i = trough_index; i < spikeSize;i++)
      {
        if(spk[i]>max)
        {
          max=spk[i];
          //    index_peak_from_baseline_post=i;
        }
      }
      peak_from_baseline_post=max-baseline;
    }
    else{
      max=trough;
      for (int i =0; i < trough_index;i++)
      {
        if(spk[i]<max)
        {
          max=spk[i];
          // index_peak_from_baseline_pre=i;
        }
      }
      peak_from_baseline_pre=max-baseline;
      max=trough;
      for (int i = trough_index; i < spikeSize;i++)
      {
        if(spk[i]<max)
        {
          max=spk[i];
          //    index_peak_from_baseline_post=i;
        }
      }
      peak_from_baseline_post=max-baseline;
    }
    
    if(sqrt(peak_from_baseline_pre*peak_from_baseline_pre)+sqrt(peak_from_baseline_post*peak_from_baseline_post)!=0){
      peak_amplitude_asymmetry=(peak_from_baseline_pre-peak_from_baseline_post)/
        (sqrt(peak_from_baseline_pre*peak_from_baseline_pre)+sqrt(peak_from_baseline_post*peak_from_baseline_post));
    }else{
      peak_amplitude_asymmetry=0;
    }
    
    
    o[0*nSpikes+i]=amplitude;
    o[1*nSpikes+i]=spike_duration;
    o[2*nSpikes+i]=first_half_duration;
    o[3*nSpikes+i]=second_half_duration;
    o[4*nSpikes+i]=peak_amplitude_asymmetry;
    
  //  Rprintf("b1:%lf b2:%lf base:%lf trou:%lf amp:%lf spikeD:%lf firstHalfD:%lf secondHalfD:%lf peakPre:%lf peakPost:%lf asym:%lf\n",
    //        b1,b2,baseline,trough,amplitude,spike_duration,first_half_duration,second_half_duration,
      //      peak_from_baseline_pre, peak_from_baseline_post,peak_amplitude_asymmetry);
    }
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
  
  // set the array to max possible size
  int* ti=(int*)malloc(size*sizeof(int));
  int num_spikes=0;
  
  int j;
  int index_keep=0;
  int smallest_trough=0;
  int latest_trough;
  for(int i = 0; i < size; i++)
  {
    // loop from a new spike or the last that was added from simultaneous spikes
    
    j=i;
    while((j+1)<size && (time[j+1]-time[i])<max_time_difference){
      //Rprintf("join with j:%d, time[j]: %d\n",j,time[j]);
      j++;
    }
    
    if(i==j){ // no simultaneous spikes to join
      if(num_spikes>0){
       if(ti[num_spikes-1]!=time[i]){ // if the spike was not already added
         ti[num_spikes]=time[i];
         latest_trough=trough[i];
         num_spikes++;
       }
      }
      else{ // this is the first spike
      ti[num_spikes]=time[i];
      latest_trough=trough[i];
      num_spikes++;
      }
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
      if(num_spikes>0){
        if(ti[num_spikes-1]+max_time_difference<time[index_keep]){ // not in the refractory of previous spike, add spike
          ti[num_spikes]=time[index_keep];
          latest_trough=trough[index_keep];
          num_spikes++;
          time[j]=time[index_keep];
          trough[j]=trough[index_keep];
          i=j-1;
        }
        else{ // in the refractory of previous spike, keep only the one with smallest trough
          if(latest_trough>trough[index_keep]) // trough is smaller for new spike, replace last added
          {
            ti[num_spikes-1]=time[index_keep];
            latest_trough=trough[index_keep];
            time[j]=time[index_keep];
            trough[j]=trough[index_keep];
            i=j-1;
          }
          else
          {
            i=j; // if we don't add the new spike and keep the already added one
          }
        }
      }
      else{ // num_spikes==0
        ti[num_spikes]=time[index_keep];
        latest_trough=trough[index_keep];
        num_spikes++;
        time[j]=time[index_keep];
        trough[j]=trough[index_keep];
        i=j-1;
        if(i<0)i=0;
      }
    }
  }
  
  SEXP out = PROTECT(allocVector(REALSXP, num_spikes));
  double* outc = REAL(out);
  for(int i = 0; i < num_spikes; i++)
    outc[i]=ti[i];
  
  free(ti);
  UNPROTECT(1);
  return(out);
}

