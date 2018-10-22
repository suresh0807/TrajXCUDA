//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################


#include "cudatools.cuh"

__global__ void SDreduceHB(float *a, float *b, int SD_store, int initHBnum)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  if(rowid < SD_store )
  {
    for(int i=0;i<initHBnum;i++)
    {
      b[rowid]+=a[i+rowid*initHBnum];
    }
    b[rowid] = b[rowid]/float(initHBnum);
  }
}


