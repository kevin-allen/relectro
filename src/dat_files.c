#include <fcntl.h> // for the open function
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include "relectro.h"

#define MAXFILENAMELENGTH 255
#define MAXNUMBERFILES 100
#define MAXBLOCKSIZE 50000000 // 50MB to work with many files at same time

// write the cwrap functions 
// get data one channels
// detect up
// detect down

SEXP group_data_file_si_get_one_channel_cwrap(SEXP file_names_r, SEXP num_channels_r, SEXP channel_no_r, SEXP start_index_r, SEXP end_index_r)
{
  // get the name of the first file
  char* data_file_names[MAXNUMBERFILES];
  int i, file_lines;
  struct group_data_file_si gdf;
  int num_samples;
  int needed_samples;

  PROTECT(file_names_r = AS_CHARACTER(file_names_r));
  file_lines = LENGTH(file_names_r);
  if(file_lines>MAXNUMBERFILES)
    {
      Rprintf("Number of files is larger than %d\n",MAXNUMBERFILES);
      UNPROTECT(1);
      return(R_NilValue);
    }

  //copy the file names into a x-element array of pointer */
  for (i=0; i<file_lines; i++) {
    if((data_file_names[i]=malloc(strlen( CHAR(STRING_ELT(file_names_r, i))) + 1))==NULL)
      {
	Rprintf("group_data_File_si_get_one_channel_cwrap(): problem allocating memory for data_file_names\n");
	UNPROTECT(1);
	return(R_NilValue);
      }
    strcpy(data_file_names[i],CHAR(STRING_ELT(file_names_r, i)));
    //    Rprintf("%s\n",data_file_names[i]);
  }
  
  if(init_group_data_file_si(&gdf,data_file_names,file_lines,INTEGER_VALUE(num_channels_r))!=0)
    {
      Rprintf("problem with init_group_data_file_si\n");
      UNPROTECT(1);
      return (R_NilValue);
    }
  num_samples=(int)gdf.num_samples_all_files;
  
  if(INTEGER_VALUE(end_index_r) > num_samples)
    {
      Rprintf("end_index_r > num_samples\n");
      UNPROTECT(1);
      return (R_NilValue);
    }

  needed_samples=INTEGER_VALUE(end_index_r)-INTEGER_VALUE(start_index_r)+1;
  SEXP out = PROTECT(allocVector(INTSXP,needed_samples));
  int* ptr = INTEGER_POINTER(out);
  
  if ((group_data_file_si_get_data_one_channel(&gdf,
  					       INTEGER_VALUE(channel_no_r),
  					       ptr,
  					       INTEGER_VALUE(start_index_r),
  					       INTEGER_VALUE(end_index_r)))!=0)
    {
      Rprintf("error reading from data file\n");
      UNPROTECT(2);
      return (R_NilValue);
    }

  // free memory
  for (i=0; i<file_lines; i++)
    free(data_file_names[i]);
  clean_group_data_file_si(&gdf);
  UNPROTECT(2);
  return(out);
}


SEXP group_data_file_si_get_group_channels_cwrap(SEXP file_names_r, SEXP num_channels_r, SEXP channels_r, SEXP num_channels_get_r, SEXP start_index_r, SEXP end_index_r)
{
  // get the name of the first file
  char* data_file_names[MAXNUMBERFILES];
  int i, file_lines;
  struct group_data_file_si gdf;
  int num_samples;
  int needed_samples;
  int* channels;
  int num_channels_get=INTEGER_VALUE(num_channels_get_r);
  if(num_channels_get<1)
  {
    Rprintf("Number of channels to get is smaller than 1\n");
    return(R_NilValue);
  }
  channels=INTEGER_POINTER(channels_r);
 
  Rprintf("num_channls_get: %d\n",num_channels_get);
  for(int i = 0; i < num_channels_get; i++)
    Rprintf("channel %d: %d\n",i,channels[i]);
  
 
 
  
  PROTECT(file_names_r = AS_CHARACTER(file_names_r));
  file_lines = LENGTH(file_names_r);
  if(file_lines>MAXNUMBERFILES)
  {
    Rprintf("Number of files is larger than %d\n",MAXNUMBERFILES);
    UNPROTECT(1);
    return(R_NilValue);
  }
  
  //copy the file names into a x-element array of pointer */
  for (i=0; i<file_lines; i++) {
    if((data_file_names[i]=malloc(strlen( CHAR(STRING_ELT(file_names_r, i))) + 1))==NULL)
    {
      Rprintf("group_data_File_si_get_group_channels_cwrap(): problem allocating memory for data_file_names\n");
      UNPROTECT(1);
      return(R_NilValue);
    }
    strcpy(data_file_names[i],CHAR(STRING_ELT(file_names_r, i)));
    //    Rprintf("%s\n",data_file_names[i]);
  }
  
  if(init_group_data_file_si(&gdf,data_file_names,file_lines,INTEGER_VALUE(num_channels_r))!=0)
  {
    Rprintf("problem with init_group_data_file_si\n");
    UNPROTECT(1);
    return (R_NilValue);
  }
  num_samples=(int)gdf.num_samples_all_files;
  
  if(INTEGER_VALUE(end_index_r) > num_samples)
  {
    Rprintf("end_index_r > num_samples\n");
    UNPROTECT(1);
    return (R_NilValue);
  }
  
  needed_samples=INTEGER_VALUE(end_index_r)-INTEGER_VALUE(start_index_r)+1;
  SEXP out = PROTECT(allocVector(INTSXP,needed_samples*num_channels_get));
  int* ptr = INTEGER_POINTER(out);
  
  Rprintf("The out array pointer is at %d\n",ptr);
  Rprintf("call to group_data_file_si_get_data_group_channels\n");
  if ((group_data_file_si_get_data_group_channels(&gdf,
                                               channels,
                                               num_channels_get,
                                               ptr, //allocated memory for the data
                                               INTEGER_VALUE(start_index_r),
                                               INTEGER_VALUE(end_index_r)))!=0)
  {
    Rprintf("error reading from data file\n");
    UNPROTECT(2);
    return (R_NilValue);
  }
  
  // free memory
  for (i=0; i<file_lines; i++)
    free(data_file_names[i]);
  clean_group_data_file_si(&gdf);
  UNPROTECT(2);
  return(out);
}

int init_group_data_file_si(struct group_data_file_si* gdf, char** file_names,int num_files,int num_channels)
{
  // function to initialize the variable and allocate memory
  int i;
  gdf->num_files=num_files;
  gdf->num_channels=num_channels;
  
  //  allocate memory for data_file_si structures
  if((gdf->file_group=(struct data_file_si*)malloc(sizeof(struct data_file_si)*gdf->num_files))==NULL)
    {
      Rprintf("init_group_data_file_si(): problem allocating memory for file_group\n");
      return 1;
    }
  // initiate all the individual files
  for (i=0; i < gdf->num_files;i++)
    {
      if((init_data_file_si(&gdf->file_group[i],file_names[i], gdf->num_channels))!=0)
  	{
  	  Rprintf("init_group_data_file_si(): problem init_data_file_si for %s\n",file_names[i]); // is this correct???????
  	  return 1;
  	}
    }
  // allocate memory for the resofs
  if((gdf->resofs=(long*)malloc(sizeof(long)*gdf->num_files))==NULL)
    {
      Rprintf("init_group_data_file_si(): problem allocating memory for resofs\n");
      return 1;
    }
  for (i=0; i < gdf->num_files;i++)
    {
      if(i==0)
  	{
  	  gdf->resofs[i]=gdf->file_group[i].num_samples_in_file;
  	}
      else
  	{
  	  gdf->resofs[i]=gdf->resofs[i-1]+gdf->file_group[i].num_samples_in_file;
  	}
    }
  // number channels all files
  gdf->num_samples_all_files=0;
  for (i=0; i < num_files;i++)
    {
      gdf->num_samples_all_files+=gdf->file_group[i].num_samples_in_file;
    }

  //Rprintf("init_group_data_file_si(): gdf->num_files: %d \n",gdf->num_files);
  //Rprintf("init_group_data_file_si(): gdf->num_channels: %d \n",gdf->num_channels);
  return 0;
}

int group_data_file_si_get_data_one_channel(struct group_data_file_si * gdf,int channel_no, int* one_channel, long int start_index, long int end_index)
{
  // function to get the data from one channel, intervals can cover more than one file
  // we do one read operation per file until we get all we need
  int* ptr;
  int file_index;
  long int num_samples_read, to_read, within_file_start_index,total_needed;
  // check that the index given make sense
  if (start_index<0)
  {
    Rprintf("group_data_file_si_get_data_one_channel(): start index is smaller than 0: %ld\n",start_index);
    return 1;
  }
  if (end_index<start_index)
  {
    Rprintf("group_data_file_si_get_data_one_channel(): end index(%ld) is smaller than start index(%ld)\n",end_index,start_index);
    return 1;
  }
  if(end_index>gdf->num_samples_all_files)
  {
    Rprintf("group_data_file_si_get_data_one_channel(): end index(%ld) is larger than the number of samples allfiles(%ld)\n",end_index,gdf->num_samples_all_files);
    return 1;
  }
  // find the file containing the beginning of data
  file_index=0;
  while(file_index+1<gdf->num_files && start_index>gdf->resofs[file_index])
  {file_index++;}
  
  
  // get the index at which to start reading operation within the selected file
  if(file_index==0)
  {
    within_file_start_index=start_index;
  }
  else
  {
    within_file_start_index=start_index-gdf->resofs[file_index-1];
  }
  
  // set variables to know when we have enough data
  num_samples_read=0;
  total_needed=end_index-start_index+1;
  while(num_samples_read<total_needed&&file_index<gdf->num_files)
  {
    to_read=total_needed-num_samples_read;
    
    // if we need more data that what is available in this file, read until the end of this file
    if(to_read > gdf->file_group[file_index].num_samples_in_file-within_file_start_index)
    {
      to_read=gdf->file_group[file_index].num_samples_in_file-within_file_start_index;
    }
    
    ptr=one_channel+num_samples_read;
    
    //  Rprintf("reading file %d from %ld to %ld\n",file_index,within_file_start_index,within_file_start_index+to_read);
    
    if((data_file_si_get_data_one_channel(&gdf->file_group[file_index], channel_no, ptr, within_file_start_index, within_file_start_index+to_read))!=0)
    {
      Rprintf("group_data_file_si_get_data_one_channel(): error reading from file %d\n",file_index);
    }
    num_samples_read+=to_read;
    within_file_start_index=0;
    file_index++;
  }
  return 0;
}

int group_data_file_si_get_data_group_channels(struct group_data_file_si* gdf,int* channels,int num_channels, 
                                               int* data, long int start_index, long int end_index)
{

  // function to get the data from several channels, intervals can cover more than one file
  // we do one read operation per file until we get all we need
  
  
  int file_index;
  long int num_samples_read, to_read, within_file_start_index,total_needed;
  // check that the index given make sense
  if (start_index<0)
  {
    Rprintf("group_data_file_si_get_data_group_channels(): start index is smaller than 0: %ld\n",start_index);
    return 1;
  }
  if (end_index<start_index)
  {
    Rprintf("group_data_file_si_get_data_group_channels(): end index(%ld) is smaller than start index(%ld)\n",end_index,start_index);
    return 1;
  }
  if(end_index>gdf->num_samples_all_files)
  {
    Rprintf("group_data_file_si_get_data_group_channels(): end index(%ld) is larger than the number of samples allfiles(%ld)\n",end_index,gdf->num_samples_all_files);
    return 1;
  }
  // find the file containing the beginning of data
  file_index=0;
  while(file_index+1<gdf->num_files && start_index>gdf->resofs[file_index])
  {file_index++;}
  
  
  // get the index at which to start reading operation within the selected file
  if(file_index==0)
  {
    within_file_start_index=start_index;
  }
  else
  {
    within_file_start_index=start_index-gdf->resofs[file_index-1];
  }
  
  
  // data is filled with all samples of one channel together, we need an array of pointer to do this
  // c1s1 c1s2 c1s3 ... c2s1 c2s2 c2s3 ...
  int** ptr;
  Rprintf("int** ptr has %d length",num_channels);
  ptr= malloc(sizeof(int*)*num_channels);
  
  
  // set variables to know when we have enough data
  num_samples_read=0;
  total_needed=end_index-start_index+1;
  while(num_samples_read<total_needed&&file_index<gdf->num_files)
  {
    to_read=total_needed-num_samples_read;
    
    // if we need more data that what is available in this file, read until the end of this file
    if(to_read > gdf->file_group[file_index].num_samples_in_file-within_file_start_index)
    {
      to_read=gdf->file_group[file_index].num_samples_in_file-within_file_start_index;
    }
    
    Rprintf("set pointers\n");
    for(int i = 0; i < num_channels; i++){
      ptr[i]=data+(total_needed*i)+(num_samples_read);   
      Rprintf("set pointer no: %d, offset: %d, pointer: %d\n",i,(total_needed*i)+(num_samples_read),ptr[i]);
    }
   
    Rprintf("num_samples_read: %d\n",num_samples_read);
    Rprintf("getting data from file %d from %ld to %ld\n",file_index,within_file_start_index,within_file_start_index+to_read);
    if((data_file_si_get_data_several_channels(&gdf->file_group[file_index], channels,num_channels, ptr, 
                                               within_file_start_index, within_file_start_index+to_read))!=0)
    {
      Rprintf("group_data_file_si_get_data_one_channel(): error reading from file %d\n",file_index);
    }
    num_samples_read+=to_read;
    within_file_start_index=0;
    file_index++;
  }

  free(ptr);
  return 0;
}



int clean_group_data_file_si(struct group_data_file_si * gdf)
{
  // function to free allocated memory
  int i;
  for(i=0; i < gdf->num_files;i++)
    {
       clean_data_file_si(&gdf->file_group[i]);
    }
  // free memory
  if(gdf->file_group!=NULL)
    free(gdf->file_group);
  if(gdf->resofs!=NULL)
    free(gdf->resofs);
  return 0;
}

int init_data_file_si(struct data_file_si *df,const char *file_name, int num_channels)
{
  /*function to initialise the data_file_si structure so that it can
   read a .dat file.
   The .dat file will remain open as long as clean_data_file_si is not called.
   This is done to speed up file operation. I will test what is the gain.
*/
  df->file_name=NULL;
  df->data_block=NULL;
  // allocate memory for df->file_name
  if((df->file_name=(char *)malloc(strlen(file_name)+1))==NULL)
    {
      Rprintf("init_data_file_si(): problem allocating memory for file_name\n");
      return 1;
    }
  strcpy(df->file_name,file_name);  // copy the file_name
  df->num_channels=num_channels; // get number of channels
  if (df->num_channels<=0)
    {
      Rprintf("init_data_file_si(): number of channels is 0 or less\n");
      return 1;
    }
  
  // give a warning if size of short is not 2 bytes on this system
  if (sizeof(short int)!=2)
    {
      Rprintf("init_data_file_si(): sizeof(short int) on this system is not 2\n");
      Rprintf("that might cause problem in reading/writing dat files\n");
    }
  
  // try to open the file
  if((df->file_descriptor=open(file_name,O_RDONLY,0))==-1)
    {
      Rprintf("init_data_file_si(): problem opening %s\n",file_name);
      return 1;
    }
  if((df->file_size=lseek(df->file_descriptor,0L,2))==-1)
    {
      Rprintf("init_data_file_si(): problem getting file size\n");
      return 1;
    }
  if(lseek(df->file_descriptor,0L,0)==-1)
    {
      Rprintf("init_data_file_si(): problem moving in the file\n");
      return 1;
    }
  if(df->file_size%(df->num_channels*sizeof(short))!=0)
    {
      Rprintf("init_data_file_si(): problem with the size of the file\n");
      Rprintf("the size should divide by num_channel*sizeof(short):%d, but is %ld\n",df->num_channels*(int)sizeof(short),(long int)df->file_size);
      return 1;
    }
  df->num_samples_in_file=df->file_size/(df->num_channels*sizeof(short));

  // calculate the size of a data_block
  df->num_samples_in_complete_block=MAXBLOCKSIZE/(df->num_channels*sizeof(short));
  df->block_size=df->num_samples_in_complete_block*(df->num_channels*sizeof(short));
  // allocate memory for data_block
  if((df->data_block=(short *)malloc(df->block_size))==NULL)
    {
      Rprintf("init_data_file_si(): problem allocating memory for data_block\n");
      return 1;
    }
  return 0;
}
int clean_data_file_si(struct data_file_si *df)
{
  /* function to free memory and close a .dat file after reading it.
   */
  // Rprintf("clean_data_file_si()\n");
  // free memory

  
  if (df->file_name!=NULL)
    free(df->file_name);
  if (df->data_block!=NULL)
    free(df->data_block);
  
  // try to close the file
  if(close(df->file_descriptor)==-1)
    {
      Rprintf("clean_data_file_si(): problem closing file descriptor %d\n",df->file_descriptor);
      return 1;
    }
  return 0;
}
int data_file_si_load_block(struct data_file_si* df, long int start_index, long int size)
{
  /* Function to read a data_block.
     Store the data in data_block member of the data_file_si structure
     start_index is in bytes from beginning of file

     assumes that the file is already open
  */
  if (size>df->block_size)
    {
      Rprintf("data_file_si_load_block(): size is larger than block_size, %ld\n",size);
      return 1;
    }
  if ( start_index < 0)
    {
      Rprintf("data_file_si_load_block(): start index < 0, %ld\n",start_index);
      return 1;
    }
  if (start_index + size > df->file_size)
    {
      Rprintf("data_file_si_load_block(): start_index+size > file_size\n");
      return 1;
    }
  if(lseek(df->file_descriptor,start_index,SEEK_SET)==-1)
    {
      Rprintf("data_file_si_load_block(): problem with lseek\n");
      return 1;
    }
  if(read(df->file_descriptor,df->data_block,size)==-1)
    {
      Rprintf("data_file_si_load_block(): problem reading the file\n");
      return 1;
    }
  return 0;
}
int data_file_si_get_data_one_channel(struct data_file_si* df, int channel_no, int* one_channel, long int start_index, long int end_index)
{
  /*function to read the data from one channel
    the index as parameters are in sample number
*/

  if(channel_no < 0)
    {
      Rprintf("data_file_si_get_data_one_channel(): channel_no < 0\n");
      return 1;
    }
  if(channel_no >= df->num_channels)
    {
      Rprintf("data_file_si_get_data_one_channel(): channel_no >= num_channels\n");
      return 1;
    }
  if (start_index<0)
    {
      Rprintf("data_file_si_get_data_one_channel(): start_index < 0\n");
      return 1;
    }
  if (start_index>df->num_samples_in_file)
    {
      Rprintf("data_file_si_get_data_one_channel(): start_index > num_samples\n");
      return 1;
    }
  if(end_index<=start_index)
    {
      Rprintf("data_file_si_get_data_one_channel(): start_index <= end_index\n");
      return 1;
    }
  if(end_index>df->num_samples_in_file)
    {
      Rprintf("data_file_si_get_data_one_channel(): end_index > num_samples\n");
      return 1;
    }
  int num_samples_to_read=end_index-start_index;
  int num_complete_blocks_to_read=num_samples_to_read/df->num_samples_in_complete_block;
  int num_blocks_to_read;
  int num_samples_incomplete_block=num_samples_to_read%df->num_samples_in_complete_block;
  long int i,j,index;
  long int start_index_bytes;
  index=0;
  if(num_samples_incomplete_block>0)
    num_blocks_to_read=num_complete_blocks_to_read+1;
  else
    num_blocks_to_read=num_complete_blocks_to_read;
  
  for (i = 0; i < num_blocks_to_read; i++)
    {
      start_index_bytes=(start_index*sizeof(short)*df->num_channels)+(df->block_size*i);
      if(i<num_complete_blocks_to_read) // complete block
	{
	  if(data_file_si_load_block(df,start_index_bytes,df->block_size)!=0)
	    {
	      Rprintf("data_file_si_get_data_one_channel(): problem loading block\n");
	      Rprintf("data_file_si_load_block(file,%ld,%d)\n",start_index_bytes,df->block_size);
	      return 1;
	    }
	  for (j = 0; j <  df->num_samples_in_complete_block; j++)
	    {
	      one_channel[index]=df->data_block[(j*df->num_channels)+(channel_no)];
	      index++;
	    }
	}
      if(i==num_complete_blocks_to_read) // smaller and last block
	{
	  if(data_file_si_load_block(df,start_index_bytes,num_samples_incomplete_block*sizeof(short)*df->num_channels))
	    {
	      Rprintf("data_file_si_get_data_one_channel(): problem loading last block\n");
	      return 1;
	    }
	  for (j = 0; j < num_samples_incomplete_block; j++)
	    {
	      one_channel[index]=df->data_block[(j*df->num_channels)+(channel_no)];
	      index++;
	    }
	}
    }
  return 0;
}

int data_file_si_get_data_several_channels(struct data_file_si* df, int* channels, int num_channels, int** ptr, long int start_index, long int end_index)
{
  if(num_channels<=0)
  {
    Rprintf("data_file_si_get_data_several_channels(): channel_no < 0\n");
    return 1;
  }
  for(int i =0; i < num_channels;i++){
    if(channels[i] < 0)
    {
      Rprintf("data_file_si_get_data_several_channels(): channel_no < 0\n");
      return 1;
    }
  }
  for(int i =0; i < num_channels;i++){
    if(channels[i] >= df->num_channels)
    {
      Rprintf("data_file_si_get_data_several_channels(): channel_no >= num_channels\n");
      return 1;
    }
  }
  if (start_index<0)
  {
    Rprintf("data_file_si_get_data_several_channels(): start_index < 0\n");
    return 1;
  }
  if (start_index>df->num_samples_in_file)
  {
    Rprintf("data_file_si_get_data_several_channels(): start_index > num_samples\n");
    return 1;
  }
  if(end_index<=start_index)
  {
    Rprintf("data_file_si_get_data_several_channels(): start_index <= end_index\n");
    return 1;
  }
  if(end_index>df->num_samples_in_file)
  {
    Rprintf("data_file_si_get_data_several_channels(): end_index > num_samples\n");
    return 1;
  }
  
  int num_samples_to_read=end_index-start_index;
  int num_complete_blocks_to_read=num_samples_to_read/df->num_samples_in_complete_block;
  int num_blocks_to_read;
  int num_samples_incomplete_block=num_samples_to_read%df->num_samples_in_complete_block;
  long int i,j,k,index;
  long int start_index_bytes;
  index=0;
  if(num_samples_incomplete_block>0)
    num_blocks_to_read=num_complete_blocks_to_read+1;
  else
    num_blocks_to_read=num_complete_blocks_to_read;
  
  Rprintf("num_blocks_to_read: %d\n",num_blocks_to_read);
  
  for(int i = 0; i < num_channels; i++)
  {
    Rprintf("ptr %d %d\n",i,ptr[i]);
  }
 
  for (i = 0; i < num_blocks_to_read; i++)
  {
    start_index_bytes=(start_index*sizeof(short)*df->num_channels)+(df->block_size*i);
    if(i<num_complete_blocks_to_read) // complete block
    {
      
      Rprintf("Read complete block %d\n",i);
      if(data_file_si_load_block(df,start_index_bytes,df->block_size)!=0)
      {
        Rprintf("data_file_si_get_data_several_channels(): problem loading block\n");
        Rprintf("data_file_si_load_block(file,%ld,%d)\n",start_index_bytes,df->block_size);
        return 1;
      }
      for (j = 0; j <  df->num_samples_in_complete_block; j++){
        for(k = 0; k < num_channels;k++)
        {
          ptr[k][index]=df->data_block[(j*df->num_channels)+(channels[k])];
        }
        index++;
      }
    }
    if(i==num_complete_blocks_to_read) // smaller and last block
    {
      Rprintf("Read incomplete block %d\n",i);
      if(data_file_si_load_block(df,start_index_bytes,num_samples_incomplete_block*sizeof(short)*df->num_channels))
      {
        Rprintf("data_file_si_get_data_several_channels(): problem loading last block\n");
        return 1;
      }
      for (j = 0; j < num_samples_incomplete_block; j++){
        for(k = 0; k < num_channels;k++)
        {
          ptr[k][index]=df->data_block[(j*df->num_channels)+(channels[k])];
        }
        index++;
      }
    }
  }
  return 0;
}  

int data_file_si_get_data_all_channels(struct data_file_si* df, short int* data, long int start_index, long int end_index)
{
  
  if (start_index<0)
    {
      Rprintf("data_file_si_get_data_all_channels(): start_index < 0\n");
      return 1;
    }
  if (start_index>df->num_samples_in_file)
    {
      Rprintf("data_file_si_get_data_all_channels(): start_index > num_samples\n");
      return 1;
    }
  if(end_index<=start_index)
    {
      Rprintf("data_file_si_get_data_all_channels(): start_index <= end_index\n");
      return 1;
    }
  if(end_index>df->num_samples_in_file)
    {
      Rprintf("data_file_si_get_data_all_channels(): end_index > num_samples\n");
      return 1;
    }
  int num_samples_to_read=end_index-start_index;
  int num_complete_blocks_to_read=num_samples_to_read/df->num_samples_in_complete_block;
  int num_blocks_to_read;
  int num_samples_incomplete_block=num_samples_to_read%df->num_samples_in_complete_block;
  int i,j,index;
  int start_index_bytes;
  index=0;
  if(num_samples_incomplete_block>0)
    num_blocks_to_read=num_complete_blocks_to_read+1;
  else
    num_blocks_to_read=num_complete_blocks_to_read;
  
  for (i = 0; i < num_blocks_to_read; i++)
    {
      start_index_bytes=(start_index*sizeof(short)*df->num_channels)+(df->block_size*i);
      if(i<num_complete_blocks_to_read) // complete block
	{
	  if(data_file_si_load_block(df,start_index_bytes,df->block_size)!=0)
	    {
	      Rprintf("data_file_si_get_data_all_channels(): problem loading block\n");
	      return 1;
	    }
	  for (j = 0; j < df->block_size/(int)sizeof(short); j++)
	    {
	      data[index]=df->data_block[j];
	      index++;
	    }
	}
      if(i==num_complete_blocks_to_read) // smaller and last block
	{
	  if(data_file_si_load_block(df,start_index_bytes,num_samples_incomplete_block*sizeof(short)*df->num_channels))
	    {
	      Rprintf("data_file_si_get_data_all_channels(): problem loading last block\n");
	      return 1;
	    }
	  for (j = 0; j < num_samples_incomplete_block*df->num_channels; j++)
	    {
	      data[index]=df->data_block[j];
	      index++;
	    }
	}
    }
  return 0;
}
