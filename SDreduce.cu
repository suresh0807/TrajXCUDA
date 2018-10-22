//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"

__global__ void SDreduce(float *a, float *b, int SD_store, int nElt1)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  if(rowid < SD_store )
  {
    for(int i=0;i<nElt1;i++)
    {
      b[rowid]+=a[i+rowid*nElt1];
    }
    b[rowid] = b[rowid]/float(nElt1);
  }
}

__global__ void SDreduce(float *a, float *b, int SD_store, int nElt1,int fairy)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  if(rowid < SD_store )
  {
    for(int i=0;i<nElt1;i++)
    {
      b[rowid]+=a[i+rowid*nElt1];
    }
  }
}

