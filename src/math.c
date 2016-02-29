#include "relectro.h"

int find_max(int num_data, int* data)
{
  /* returns the maximum value in an array */
  int max=data[0];
  for (int i =1; i < num_data; i++)
    {
      if (data[i]>max)
	{
	  max=data[i];
	}
    }
  return max;
}
