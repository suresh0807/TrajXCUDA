//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################





#include "cudatools.cuh"

__global__ void Xdist_calc(float *A,int *density,int nstruct,int natoms,int zsplit, float *ztick, int Dirn)
{
  int rowid=threadIdx.x + blockIdx.x * blockDim.x;
  int colid=threadIdx.y + blockIdx.y * blockDim.y;
  if(rowid < nstruct && colid <natoms)
  {
        for(int m=0;m<zsplit;m++)
        {
        if( A[colid*3+(rowid*natoms*3)+Dirn] >= ztick[m*2] && A[colid*3+(rowid*natoms*3)+Dirn] < ztick[m*2+1])
          {
            atomicAdd(&density[m],1);
            break;
  //          density[m+l*zsplit+k*ysplit*zsplit]+=1;
          }
        }
   }
}


__global__ void Xdist_calc(float *A, float *B,float *tetradensity, int *density,int nstruct,int natoms,int zsplit, float *ztick, int Dirn)
{
  int rowid=threadIdx.x + blockIdx.x * blockDim.x;
  int colid=threadIdx.y + blockIdx.y * blockDim.y;
  if(rowid < nstruct && colid <natoms)
  {
        for(int m=0;m<zsplit;m++)
        {
        if( A[colid*3+(rowid*natoms*3)+Dirn] >= ztick[m*2] && A[colid*3+(rowid*natoms*3)+Dirn] < ztick[m*2+1])
          {
	    atomicAdd(&density[m],1);
            atomicAdd(&tetradensity[m],B[colid+rowid*natoms]);
            break;
  //          density[m+l*zsplit+k*ysplit*zsplit]+=1;
          }
        }
   }
}