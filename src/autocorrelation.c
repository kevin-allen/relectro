#include "relectro.h"

// function from libelectro, remove all c++ components: 
// 3 assert statements removed and histo is now int *
void autocorrelation_one_cell(int clu_no, // cell of interest
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
  int time_diff;
  int index,k;
  int min = 0 - window_size/2;
  int max = 0 + window_size/2;
  double interval = (double)(max - min)/histo_size; // bin_size
  
  // set to 0
  for (int i = 0 ; i < histo_size; i++)
    {
      histo[i]=0;
    }
   /// for every interval
  for(int inter = 0; inter < interval_lines; inter++) 
    {
      for (int j = start_interval_index[inter]; j <= end_interval_index[inter]&&j<res_lines; j++){ // for every spikes within interval
	
	  // if cell of interest fires
	  if (clu[j] == clu_no)
	    {
	      // check if there is a spike from the same cell in the window_size
	      // check backwark
	      for (k = j-1; (k >= start_interval_index[inter]&&k<res_lines) && (res[k] > res[j] - window_size/2); k--)
		{
		  if (clu[k] == clu_no)
		    {
		      time_diff = res[k] - res[j];
		      index = (int)(time_diff/interval + histo_size/2);
		      histo[index]++;
		    }
		}
	      // check forward
	      for (k = j+1; k <= end_interval_index[inter] && k < res_lines  && (res[k] < res[j] + window_size/2); k++)
		{
		  if (clu[k] == clu_no)
		    {
		      time_diff = res[k] - res[j];
		      index =(int)( time_diff/interval + histo_size/2);
		      histo[index]++;
		    }
		}
	    }
	}
    }
  return;
 }

void autocorrelation_one_cell_probability(int clu_no, // cell of interest
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
  int time_diff;
  int index,k;
  int min = 0 - window_size/2;
  int max = 0 + window_size/2;
  double interval = (double)(max - min)/histo_size; // bin_size
  int count=0;
  // set to 0
  for (int i = 0 ; i < histo_size; i++)
    {
      histo[i]=0;
    }
   /// for every interval
  for(int inter = 0; inter < interval_lines; inter++) 
    {
      for (int j = start_interval_index[inter]; j <= end_interval_index[inter]&& j< res_lines; j++) // for every spikes within interval
	{
	  // if cell of interest fires
	  if (clu[j] == clu_no)
	    {
	      count++;
	      // check if there is a spike from the same cell in the window_size
	      // check backwark
	      for (k = j-1; (k >= start_interval_index[inter] && k < res_lines) && (res[k] > res[j] - window_size/2); k--)
		{
		  if (clu[k] == clu_no)
		    {
		      time_diff = res[k] - res[j];
		      index = (int)(time_diff/interval + histo_size/2);
		      histo[index]++;
		    }
		}
	      // check forward
	      for (k = j+1; k <= end_interval_index[inter] && k < res_lines && (res[k] < res[j] + window_size/2); k++)
		{
		  if (clu[k] == clu_no)
		    {
		      time_diff = res[k] - res[j];
		      index =(int)( time_diff/interval + histo_size/2);
		      histo[index]++;
		    }
		}
	    }
	}
    }
  if(count!=0){
  for(int i = 0; i < histo_size;i++)
    histo[i]=histo[i]/count;
  }
  return;
 }

// function that loops for all our cells
void autocorrelation(int* cell_list, // cells of interest
		    int cell_list_lines,
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
  for(int i = 0; i < cell_list_lines;i++)
    {
      ptr=histo+(histo_size*i); // pointer to a single histogram
      autocorrelation_one_cell(cell_list[i],
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
void autocorrelation_probability(int* cell_list, // cells of interest
				 int cell_list_lines,
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
  for(int i = 0; i < cell_list_lines;i++)
    {
      ptr=histo+(histo_size*i); // pointer to a single histogram
      autocorrelation_one_cell_probability(cell_list[i],
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

// wrapper to call autocorrelation from .Call()
// same name as c function with _cwrap suffix
// all object types are SEXP
// argument names have the suffix _r
// should deal with count and probability
SEXP autocorrelation_cwrap(SEXP cell_list_r, 
			   SEXP cell_list_lines_r,
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
  int probability;
  
  // protect R object created in c code so that R does not delete them
  PROTECT(cell_list_r=AS_INTEGER(cell_list_r));
  PROTECT(clu_r=AS_INTEGER(clu_r));
  PROTECT(res_r=AS_INTEGER(res_r));
  PROTECT(start_interval_index_r=AS_INTEGER(start_interval_index_r));
  PROTECT(end_interval_index_r=AS_INTEGER(end_interval_index_r));
 
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

  // might want to check the arguments here
  if(probability!=0&probability!=1)
    {
      Rprintf("probability needs to be 0 or 1 but was %d\n",probability);
      UNPROTECT(5);
      return(R_NilValue);
    }
  if(probability==0)
    {
      // prepare a SEXP object with the data from the histogram
      int all_histo_length=cell_list_lines*histo_size;
      SEXP out = PROTECT(allocVector(INTSXP, all_histo_length));
      // allocate memory for histo
      int* histo = (int*)malloc(histo_size*cell_list_lines*sizeof(int));
      // call the c function from electrophys library
      autocorrelation(cell_list, // cells of interest
  		      cell_list_lines,
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
      UNPROTECT(6);
      // return the histogram
      return(out);
    }
  if(probability==1)
    {
      int all_histo_length=cell_list_lines*histo_size;
      SEXP out = PROTECT(allocVector(REALSXP, all_histo_length));
      // allocate memory for histo
      double* histo = (double*)malloc(histo_size*cell_list_lines*sizeof(double));

      autocorrelation_probability(cell_list, // cells of interest
  		      cell_list_lines,
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
      UNPROTECT(6);
      return(out);
    }
  // unprotect memory
  UNPROTECT(5);
  //
  return(R_NilValue);
}
