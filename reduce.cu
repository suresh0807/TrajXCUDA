//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"

__global__ void reduce(float *a, float *b, int bin, int nstruct)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  if(rowid < bin )
  {
    for(int i=0;i<nstruct;i++)
    {
      b[rowid]+=a[rowid+(i*bin)];
    }
  }
}

