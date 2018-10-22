//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"

__global__ void reducef(float *a, float *b, int bin, int nstruct, int nElt1)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int colid=threadIdx.y + blockIdx.y*blockDim.y;
  if(rowid < bin && colid < nstruct)
  {
    for(int i=0;i<nElt1;i++)
    {
      b[rowid+colid*bin]+=a[rowid+(i*bin)+(colid*nElt1*bin)];
    }
  }
}

