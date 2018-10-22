//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"

float maximum(float *a,int natoms,int stride)
{
  int j;
  float max;
  max=a[stride];
  for(j=1;j<natoms;j++)
  {
    if(a[j*3+stride] > max)
    {
      max=a[j*3+stride];
    }
  }
return max;
}

