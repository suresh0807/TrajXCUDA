//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################



#include "cudatools.cuh"

__global__ void VAF_calc(float *A,float *B,float *C,int *D,int SD_store,int nElt1, int i, int j,int origins,int skips,int whichwater)
{
  int rowid=threadIdx.x + blockIdx.x * blockDim.x;
  //int colid=threadIdx.y + blockIdx.y * blockDim.y;
  int widid=threadIdx.z + blockIdx.z * blockDim.z;
  float chk,z;
  if(rowid < nElt1  &&  widid < origins/skips)
  {

    z=0.0;
    for(int k=0; k< 3 ; k++)
    {
      chk=A[rowid*3+(widid*skips*nElt1*3)+(i*nElt1*3*SD_store)+k]*B[rowid*3+(j*nElt1*3)+(widid*skips*nElt1*3)+(i*nElt1*3*SD_store)+k];
      z+=chk;
    }
//     printf("%f \n",C[rowid+(j*nElt1)+(widid*nElt1*SD_store)]);
    //if(D[rowid+(widid*skips*nElt1)+(i*nElt1*SD_store)] ==whichwater && D[rowid+(j*nElt1)+(widid*skips*nElt1)+(i*nElt1*SD_store)] ==whichwater)
    if(D[rowid+(widid*skips*nElt1)+(i*nElt1*SD_store)] ==whichwater)
    {
    C[rowid+(j*nElt1)+widid*nElt1*SD_store] += z ;
    }
    else if(whichwater==2)
    {
     C[rowid+(j*nElt1)+widid*nElt1*SD_store] += z ; 
    }
  }
}

