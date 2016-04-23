#include "relectro.h"

void set_res_outside_interval_to_minus_one(int interval_lines,
					   int* start_interval_index,
					   int* end_interval_index,
					   int res_lines,
					   int* res)
{
  /*set the value in the res array to -1 if outside the intervals */
  // values at start_interval_index and end_interval_index are within intervals
  for(int i =0; i < interval_lines; i++)
    {
      if (i==0)// first interval, remove what was before the start
	{
	  for (int j = 0; j < start_interval_index[i]; j++)
	    {
	      res[j]=-1;
	    }
	}
      else // not first interval, remove from end of previous to start of current interval
	{
	  // end_interval_index[i-1]+1, is the first res value after the previous interval
	  for (int j = end_interval_index[i-1]+1; j < start_interval_index[i]; j++)
	    {
	      res[j]=-1;
	    }
	}
      if(i==interval_lines-1) // last interval, remove all after
	{
	  for (int j = end_interval_index[i]+1; j < res_lines; j++)
	    {
	      res[j]=-1;
	    }
	}
    }
  return;
}


// function from libelectro, remove all c++ components: 
// 3 assert statements removed and histo is now int *
void crosscorrelation_one_cell(int clu_no1, // cell of interest1
			       int clu_no2, // cell of interest2
			       int* clu,  
			       int* res, 
			       int res_lines,
			       int* histo, // pointer to one histogram
			       int histo_size, // histo_size
			       int window_size, // size of the window_size in res value
			       int* start_interval_index,
			       int* end_interval_index,
			       int interval_lines)
{
  /* make a crosscorrelation for a pair of cells, within intervals only */
  int time_diff;
  int index,k;
  int min = 0 - window_size/2;
  int max = 0 + window_size/2;
  double interval = (double)(max - min)/histo_size; // bin_size
  for (int i = 0; i < histo_size; i++)
    {
      histo[i] = 0;
    }
  for(int inter = 0; inter < interval_lines; inter++)  // for every interval
    { 
      for (int j = start_interval_index[inter]; j < end_interval_index[inter]; j++) // for every spikes
	{
	  // if cell 1 fires
	  if (clu[j] == clu_no1)
	    {
	      // check if there is a spike from cell2 in the window
	      // check backwark, within the interval and within the time window
	      for (k = j-1; (k >= start_interval_index[inter]) && (res[k] > res[j] - window_size/2); k--)
		{
		  if (clu[k] == clu_no2)
		    {
		      time_diff = res[k] - res[j];
		      index = (int)(time_diff/interval + histo_size/2);
		      histo[index]++;
		    }
		}
	      // check forward within the interval and within the time window
	      for (k = j+1; k < end_interval_index[inter] && (res[k] < res[j] + window_size/2); k++)
		{
		  if (clu[k] == clu_no2)
		    {
		      time_diff = res[k] - res[j];
		      index = (int)(time_diff/interval + histo_size/2);
		      histo[index]++;
		    }
		  
		}
	    }
	}
    }
  
}

void crosscorrelation_one_cell_probability(int clu_no1, // cell of interest1
					   int clu_no2, // cell of interest2
					   int* clu,  
					   int* res, 
					   int res_lines,
					   double* histo, // pointer to one histogram
					   int histo_size, // histo_size
					   int window_size, // size of the window_size in res value
					   int* start_interval_index,
					   int* end_interval_index,
					   int interval_lines)
{
  /* calculate the probability of crosscorrelation events between pair of cells within intervals */
  int count = 0;
  int time_diff;
  int index,k;
  int min = 0 - window_size/2;
  int max = 0 + window_size/2;
  double interval = (double)(max - min)/histo_size; // bin_size
  for (int i = 0; i < histo_size; i++)
    {
      histo[i] = 0;
    }
  for(int inter = 0; inter < interval_lines; inter++)  // for every interval
    { 
      for (int j = start_interval_index[inter]; j < end_interval_index[inter]; j++) // for every spikes
	{
	  // if cell 1 fires
	  if (clu[j] == clu_no1)
	    {
	      count++;
	      // check if there is a spike from cell2 in the window
	      // check backwark, within the interval and within the time window
	      for (k = j-1; (k >= start_interval_index[inter]) && (res[k] > res[j] - window_size/2); k--)
		{
		  if (clu[k] == clu_no2)
		    {
		      time_diff = res[k] - res[j];
		      index = (int)(time_diff/interval + histo_size/2);
		      histo[index]++;
		    }
		}
	      
	      // check forward within the interval and within the time window
	      for (k = j+1; k < end_interval_index[inter] && (res[k] < res[j] + window_size/2); k++)
		{
		  if (clu[k] == clu_no2)
		    {
		      time_diff = res[k] - res[j];
		      index = (int)(time_diff/interval + histo_size/2);
		      histo[index]++;
		    }
		}
	    }
	}
    }
  if(count>0)
    {
      for (int i = 0; i < histo_size; i++)
	{
	  histo[i] = histo[i]/count;
	}
    }
  else
    {
      for (int i = 0; i < histo_size; i++)
	{
	  histo[i] = 0;
	}
    }
}

// function that loops for all our cells
void crosscorrelation(int* clu1, // cells of interest
		      int* clu2,
		      int cell_pair_lines,
		      int* clu,  
		      int* res, 
		      int res_lines,
		      int* histo, // pointer to one histogram
		      int histo_size, // histo_size
		      int window_size, // size of the window_size in res value
		      int* start_interval_index,
		      int* end_interval_index,
		      int interval_lines)
{

  int* ptr;
  for(int i = 0; i < cell_pair_lines;i++)
    {
      ptr=histo+(histo_size*i); // pointer to a single histogram
      crosscorrelation_one_cell(clu1[i],
				clu2[i],
				clu,
				res,
				res_lines,
				ptr, // pointer to one histogram
				histo_size, // histo_size
				window_size, // size of the window_size in res value
				start_interval_index,
				end_interval_index,
				interval_lines);
    }
  return;
}


// function that loops for all our cells
void crosscorrelation_probability(int* clu1, // cells of interest
				  int* clu2,
				  int cell_pair_lines,
				  int* clu,  
				  int* res, 
				  int res_lines,
				  double* histo, // pointer to one histogram
				  int histo_size, // histo_size
				  int window_size, // size of the window_size in res value
				  int* start_interval_index,
				  int* end_interval_index,
				  int interval_lines)
{
  double* ptr;
  for(int i = 0; i < cell_pair_lines;i++)
    {
      ptr=histo+(histo_size*i); // pointer to a single histogram
      crosscorrelation_one_cell_probability(clu1[i],
					    clu2[i],
					    clu,
					    res,
					    res_lines,
					    ptr, // pointer to one histogram
					    histo_size, // histo_size
					    window_size, // size of the window_size in res value
					    start_interval_index,
					    end_interval_index,
					    interval_lines);
    }
  return;
}

// wrapper to call crosscorrelation from .Call()
// same name as c function with _cwrap suffix
// all object types are SEXP
// argument names have the suffix _r
// should deal with count and probability
SEXP crosscorrelation_cwrap(SEXP clu1_r, 
			    SEXP clu2_r,
			    SEXP cell_pair_lines_r,
			    SEXP clu_r,
			    SEXP res_r,
			    SEXP res_lines_r,
			    SEXP histo_size_r, // histo_size in bins
			    SEXP window_size_r, // size of the window_size in res value
			    SEXP start_interval_index_r,
			    SEXP end_interval_index_r,
			    SEXP interval_lines_r,
			    SEXP probability_r) // flag: 0 = count, 1 = probability
{ 
  // transform the SEXP object in their correct c types
  // create the list of variable of correct c types
  int* clu1;
  int* clu2;
  int cell_pair_lines;
  int* clu;
  int* res;
  int res_lines;
  int histo_size; // histo_size
  int window_size; // size of the window_size in res value
  int* start_interval_index;
  int* end_interval_index;
  int interval_lines;
  int probability;
  
  // protect R object created in c code so that R does not delete them
  PROTECT(clu1_r=AS_INTEGER(clu1_r));
  PROTECT(clu2_r=AS_INTEGER(clu2_r));
  PROTECT(clu_r=AS_INTEGER(clu_r));
  PROTECT(res_r=AS_INTEGER(res_r));
  PROTECT(start_interval_index_r=AS_INTEGER(start_interval_index_r));
  PROTECT(end_interval_index_r=AS_INTEGER(end_interval_index_r));
 
  /* // coersion */
  clu1=INTEGER_POINTER(clu1_r);
  clu2=INTEGER_POINTER(clu2_r);
  cell_pair_lines=INTEGER_VALUE(cell_pair_lines_r); //get an integer
  clu=INTEGER_POINTER(clu_r);
  res=INTEGER_POINTER(res_r);
  res_lines=INTEGER_VALUE(res_lines_r);
  histo_size=INTEGER_VALUE(histo_size_r);
  window_size=INTEGER_VALUE(window_size_r);
  start_interval_index=INTEGER_POINTER(start_interval_index_r);
  end_interval_index=INTEGER_POINTER(end_interval_index_r);
  interval_lines=INTEGER_VALUE(interval_lines_r);
  probability=INTEGER_VALUE(probability_r);


  // might want to check the arguments here
  if(probability!=0&probability!=1)
    {
      Rprintf("probability needs to be 0 or 1 but was %d\n",probability);
      UNPROTECT(6);
      return(R_NilValue);
    }
  if(probability==0)
    {
      // prepare a SEXP object with the data from the histogram
      int all_histo_length=cell_pair_lines*histo_size;
      SEXP out = PROTECT(allocVector(INTSXP, all_histo_length));
      // allocate memory for histo
      int* histo = (int*)malloc(all_histo_length*sizeof(int));
      // call the c function from electrophys library
      crosscorrelation(clu1,
		       clu2,
		       cell_pair_lines,
		       clu,
		       res,
		       res_lines,
		       histo, // pointer to all  histograms
		       histo_size, // histo_size
		       window_size, // size of the window_size in res value
		       start_interval_index,
		       end_interval_index,
		       interval_lines);
      // copy the results in a SEXP
      for(int i = 0; i < all_histo_length;i++)
  	INTEGER(out)[i]=histo[i]; 
      //free memory for histo
      free(histo);
      UNPROTECT(7);
      // return the histogram
      return(out);
    }
  if(probability==1)
    {
      int all_histo_length=cell_pair_lines*histo_size;
      SEXP out = PROTECT(allocVector(REALSXP, all_histo_length));
      // allocate memory for histo
      double* histo = (double*)malloc(all_histo_length*sizeof(double));
      crosscorrelation_probability(clu1,
				   clu2,
				   cell_pair_lines, // cells of interest
				   clu,
				   res,
				   res_lines,
				   histo, // pointer to all  histograms
				   histo_size, // histo_size
				   window_size, // size of the window_size in res value
				   start_interval_index,
				   end_interval_index,
				   interval_lines);

      for(int i = 0; i < all_histo_length;i++)
  	REAL(out)[i]=histo[i];
      free(histo);
      UNPROTECT(7);
      return(out);
    }
  // unprotect memory
  UNPROTECT(6);
  //
  return(R_NilValue);
}



SEXP crosscorrelationEvents_cwrap(SEXP cell_list_r, 
				   SEXP cell_list_lines_r,
				   SEXP clu_r,
				   SEXP res_r,
				   SEXP res_lines_r,
				   SEXP histo_size_r, // histo_size in bins
				   SEXP window_size_r, // size of the window_size in res value
				   SEXP start_interval_index_r,
				   SEXP end_interval_index_r,
				   SEXP interval_lines_r,
				   SEXP probability_r, // flag: 0 = count, 1 = probability
				   SEXP events_r,
				   SEXP events_lines_r) 
{
 
  // transform the SEXP object in their correct c types
  // create the list of variable of correct c types
  int* cell_list;
  int cell_list_lines;
  int* clu;
  int* res;
  int res_lines;
  int histo_size; // histo_size
  int window_size; // size of the window_size in res value
  int* start_interval_index;
  int* end_interval_index;
  int interval_lines;
  int* events;
  int events_lines;
  int probability;
  
  // protect R object created in c code so that R does not delete them
  PROTECT(cell_list_r=AS_INTEGER(cell_list_r));
  PROTECT(clu_r=AS_INTEGER(clu_r));
  PROTECT(res_r=AS_INTEGER(res_r));
  PROTECT(start_interval_index_r=AS_INTEGER(start_interval_index_r));
  PROTECT(end_interval_index_r=AS_INTEGER(end_interval_index_r));
  PROTECT(end_interval_index_r=AS_INTEGER(end_interval_index_r));
  PROTECT(events_r=AS_INTEGER(events_r));
 
  /* // coersion */
  cell_list=INTEGER_POINTER(cell_list_r);
  cell_list_lines=INTEGER_VALUE(cell_list_lines_r); //get an integer
  clu=INTEGER_POINTER(clu_r);
  res=INTEGER_POINTER(res_r);
  res_lines=INTEGER_VALUE(res_lines_r);
  histo_size=INTEGER_VALUE(histo_size_r);
  window_size=INTEGER_VALUE(window_size_r);
  start_interval_index=INTEGER_POINTER(start_interval_index_r);
  end_interval_index=INTEGER_POINTER(end_interval_index_r);
  interval_lines=INTEGER_VALUE(interval_lines_r);
  probability=INTEGER_VALUE(probability_r);
  events=INTEGER_POINTER(events_r);
  events_lines=INTEGER_VALUE(events_lines_r);

  // might want to check the arguments here
  if(probability!=0&probability!=1)
    {
      Rprintf("probability needs to be 0 or 1 but was %d\n",probability);
      UNPROTECT(7);
      return(R_NilValue);
    }

  if(events_lines<1)
    {
      Rprintf("events_lines:%d\n",events_lines);
      UNPROTECT(7);
      return(R_NilValue);
    }


  // prepare a SEXP object with the data from the histogram
  int all_histo_length=cell_list_lines*histo_size;
  SEXP out = PROTECT(allocVector(REALSXP, all_histo_length));
  // allocate memory for histo
  double* histo = (double*)malloc(histo_size*cell_list_lines*sizeof(double));

  cross_correlation_events(cell_list,
			   cell_list_lines,
			   clu,  
			   res, 
			   res_lines,
			   events,
			   events_lines,
			   histo, // pointer to all histograms
			   histo_size, // one histo size
			   window_size, // size of the window_size in res value
			   start_interval_index, 
			   end_interval_index, 
			   interval_lines);



  if(probability==1)
    for(int i = 0; i < all_histo_length;i++)
      histo[i]/events_lines;


  // copy the results in a SEXP
  for(int i = 0; i < all_histo_length;i++)
    REAL(out)[i]=histo[i]; 
    
  
  //free memory for histo
  free(histo);
  UNPROTECT(8);
  // return the histogram
  return(out);
  
}




void cross_correlation_events(int* cells,
			      int cell_lines,
			      int* clu,  
			      int* res, 
			      int res_lines,
			      int* events,
			      int events_lines,
			      double* histo, // pointer to all histograms
			      int histo_size, // one histo size
			      int window_size, // size of the window_size in res value
			      int* start_interval_index, 
			      int* end_interval_index, 
			      int interval_lines) 
{
  /* make crosscorrelation histograms for a list of cells, spikes around list of events */
  int* start_events_index;
  int* end_events_index;
  int* start_events;
  int* end_events;
  double* one_histo;
  int* res_mod; // modified
  int time_diff;
  int index;
  int min = 0 - window_size/2;
  int max = 0 + window_size/2;
  double interval = (max-min)/(double)histo_size; // bin_size
  // get the res_index for beginning and end of the events windows

  
  start_events =(int*)malloc(events_lines*sizeof(int));
  end_events = (int*)malloc(events_lines*sizeof(int));
  start_events_index = (int*)malloc(events_lines*sizeof(int));
  end_events_index = (int*)malloc(events_lines*sizeof(int));
  res_mod = (int*)malloc(res_lines*sizeof(int));

  // copy res into res_mod
  for (int i = 0; i < res_lines; i++)
    {
      res_mod[i]=res[i];
    }

  // set the res value outside the interval_index to -1 
  // done so we dont care about these spikes outside intervals
  set_res_outside_interval_to_minus_one(interval_lines, start_interval_index,end_interval_index,res_lines, res_mod);
  
  // calculate the beginning and end of event windows
  for (int i = 0; i < events_lines; i++)
    {
      start_events[i]=events[i]-window_size/2;
      end_events[i]=events[i]+window_size/2;
    }
  
  for(int i = 0; i < events_lines;i++)
    {
      if(start_events[i]<0)
	start_events[i]=0;
      // the end will be check in res_index_for_intervals
    }

  // get the index for the event windows
  res_index_for_intervals(&events_lines,start_events,end_events,res_lines,res,start_events_index,end_events_index);
  
  // set histo to 0
  for (int i = 0; i < cell_lines; i++)
    {
      one_histo=histo+(histo_size*i);
      for (int j = 0 ; j < histo_size; j++)
	{
	  one_histo[j]=0;
	}
    }

  // Rprintf("Number of events: %d\n",events_lines);
  //Rprintf("Number of cells: %d\n",cell_lines);

  /// for every events
  for(int inter = 0; inter < events_lines; inter++) 
    {
      //  Rprintf("Events[inter]:%d, Inter: %d\n",events[inter],inter);
      
      for (int i = 0; i < cell_lines; i++)
	{
	  one_histo=histo+(histo_size*i);

	  for (int j = start_events_index[inter]; j <= end_events_index[inter]; j++) // for every spikes within events
	    {
	      // if cell of interest fires, and the res value is within intervals of interest (not -1)
	      if (clu[j] == cells[i]&& res_mod[j]!=-1)
		{
		  time_diff = res[j] - events[inter];
		  index = (int)(time_diff/interval+histo_size/2);
		  //	  Rprintf("res:%d, events:%d, time_diff:%d, index:%d\n",res[j],events[inter],time_diff,index);
		  one_histo[index]++;
		}
	    }
	}
    }
  free(start_events_index);
  free(end_events_index);
  free(start_events);
  free(end_events);
  free(res_mod);
  return ;
}
