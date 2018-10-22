//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"

__global__ void SDreducefHB(float *a, float *b, int initHBnum, int SD_store, int origins, int skips)
{
  int rowid=threadIdx.x + blockIdx.x * blockDim.x;
  int colid=threadIdx.y + blockIdx.y * blockDim.y;
 // int widid=threadIdx.z + blockIdx.z * blockDim.z;
  if(rowid < initHBnum && colid < SD_store)
  {
    for(int i = 0; i < origins/skips; i++)
    {
      //printf("%f \n", a[rowid+(colid*nElt1)+(i*nElt1*SD_store)]); 
      b[rowid+colid*initHBnum] += a[rowid+(colid*initHBnum)+(i*initHBnum*SD_store)];
    }
    b[rowid+colid*initHBnum] = b[rowid+colid*initHBnum]/(float(origins)/float(skips));
  }
}

