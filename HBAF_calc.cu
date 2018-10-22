//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################



#include "cudatools.cuh"

__global__ void HBAF_calc(float *A,float *B,float *C,int SD_store,int initHBnum, int i, int j,int origins,int skips)
{
  int rowid=threadIdx.x + blockIdx.x * blockDim.x;
  //int colid=threadIdx.y + blockIdx.y * blockDim.y;
  int widid=threadIdx.z + blockIdx.z * blockDim.z;
  float z;
  if(rowid < initHBnum   && widid < SD_store)
  {

      z=A[rowid]*B[rowid+(widid*initHBnum)];
      
//     printf("%f \n",C[rowid+(j*nElt1)+(widid*nElt1*SD_store)]);
//    if(D[rowid+(widid*skips*nElt1)+(i*nElt1*SD_store)] ==whichwater && D[rowid+(j*nElt1)+(widid*skips*nElt1)+(i*nElt1*SD_store)] ==whichwater)
//    {
      C[rowid+widid*initHBnum] += z ;
      //atomicAdd(C[rowid+widid*initHBnum],z);
//    }
//    else if(whichwater==2)
 //   {
//     C[rowid+(j*nElt1)+widid*nElt1*SD_store] += z ; 
//    }
  }
}
