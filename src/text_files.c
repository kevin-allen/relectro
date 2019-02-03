#include "relectro.h"
#include<stdio.h>
#include<stdlib.h>

int file_lines(const char* file_name)
{
  // count the number of lines in the file called filename                                    
  FILE *fp = fopen(file_name,"r");
  int ch=0;
  int lines=0;
  if (fp == NULL){
    Rprintf("pf ==NULL\n");
    return 0;
  }
  while(!feof(fp))
  {
    ch = fgetc(fp);
    if(ch == '\n')
    {
      lines++;
    }
  }
  fclose(fp);
  return lines;
}
int read_one_column_int_file(const char* file_name,int* data, int lines)
{
  FILE *fp = fopen(file_name,"r");
  if (!fp) {
    Rprintf("problem opening the file\n");
    return 1;
  }
  for(int i = 0; i < lines; i++)
  {
    if(fscanf(fp,"%d\n",&data[i])<1){
      Rprintf("problem reading the file\n");
      return 1;
    }
  }
  fclose(fp);
  
  FILE *fp = fopen(file_name,"r");
  if (!fp) {
    Rprintf("problem opening %s\n",file_name);
    return (R_NilValue);
  }
  
  // first line is number of columns in subsequent lines
  if(fscanf(fp,"%d\n",&ncol)<1){
    printf("problem reading %s\n",file_name);
    return (R_NilValue);
  }
  if(ncol<1){
    Rprintf("read_fet_file_cwrap, ncol is < 1, %s\n",file_name);
    return (R_NilValue);
  }
  lines=lines-1;





  
  return 0;
}
SEXP read_one_column_int_file_cwrap(SEXP file_name_r)
{
  const char* file_name = CHAR(asChar(file_name_r));
  int lines;
  lines=file_lines(file_name);
  SEXP out = PROTECT(allocVector(INTSXP,lines));
  int* ptr = INTEGER_POINTER(out);
  if(read_one_column_int_file(file_name,ptr,lines)!=0){
    Rprintf("problem reading %s\n",file_name);
    UNPROTECT(1);
    return(R_NilValue);
  }
  UNPROTECT(1);
  return(out);
}
SEXP read_fet_file_cwrap(SEXP file_name_r)
{
  const char* file_name = CHAR(asChar(file_name_r));
  int lines;
  lines=file_lines(file_name);
  int ncol;
  FILE *fp = fopen(file_name,"r");
  if (!fp) {
    Rprintf("problem opening %s\n",file_name);
    return (R_NilValue);
  }
  
 // first line is number of columns in subsequent lines
  if(fscanf(fp,"%d\n",&ncol)<1){
    Rprintf("problem reading %s\n",file_name);
    return (R_NilValue);
  }
  if(ncol<1){
    Rprintf("read_fet_file_cwrap, ncol is < 1, %s\n",file_name);
    return (R_NilValue);
  }
  lines=lines-1;
  SEXP out = PROTECT(allocMatrix(INTSXP,lines,ncol));
  int* ptr = INTEGER_POINTER(out);
  for(int i = 0; i < lines; i++)
    for(int j = 0; j < ncol; j++)
      if(fscanf(fp,"%d",&ptr[j*lines+i])<1){
        Rprintf("problem reading %d element of %s\n",i,file_name);
        UNPROTECT(1);
        return (R_NilValue);
      }
  fclose(fp);
  UNPROTECT(1);
  return(out);
}
