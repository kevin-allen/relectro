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
    printf("pf ==NULL\n");
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
    printf("problem opening the file\n");
    return 1;
  }
  int ret;
  for(int i = 0; i < lines; i++)
  {
    if(fscanf(fp,"%d\n",&data[i])<1){
      printf("problem reading the file\n");
      return 1;
    }
  }
  fclose(fp);
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
    printf("problem reading %s\n",file_name);
    UNPROTECT(1);
    return(R_NilValue);
  }
  UNPROTECT(1);
  return(out);
}