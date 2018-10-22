//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"

float minimum(float *a,int natoms,int stride)
{
  int j;
  float min;
  min=a[stride];
  for(j=1;j<natoms;j++)
  {
    if(a[j*3+stride] < min)
    {
      min=a[j*3+stride];
    }
  }
return min;
}

