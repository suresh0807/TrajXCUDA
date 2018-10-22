//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




//the matrix fills Z-(vacuum) axis first and then Y and then x at last.

#include "cudatools.cuh"

__global__ void density_calc(float *A,int *density,int nstruct,int natoms,int xsplit,int ysplit, int zsplit, float *xtick,float *ytick, float *ztick)
{
  int rowid=threadIdx.x + blockIdx.x * blockDim.x;
  int colid=threadIdx.y + blockIdx.y * blockDim.y;
  if(rowid < nstruct && colid <natoms)
  {
    for(int k=0;k<xsplit;k++)
    {
      for(int l=0;l<ysplit;l++)
      {
        for(int m=0;m<zsplit;m++)
        {
        if(A[colid*3+rowid*natoms*3] >= xtick[k*2] && A[colid*3+rowid*natoms*3] < xtick[k*2+1] && A[colid*3+(rowid*natoms*3)+1] >= ytick[l*2] && A[colid*3+(rowid*natoms*3)+1] < ytick[l*2+1] && A[colid*3+(rowid*natoms*3)+2] >= ztick[m*2] && A[colid*3+(rowid*natoms*3)+2] < ztick[m*2+1])
          {
            atomicAdd(&density[m+l*zsplit+k*ysplit*zsplit],1);
	    break;
          }
           
        }
         
      }
      
    }
    }
}

__global__ void density_calc(float *A,int *density,int *exch,int nstruct,int natoms,int xsplit,int ysplit, int zsplit, float *xtick,float *ytick, float *ztick)
{
  int rowid=threadIdx.x + blockIdx.x * blockDim.x;
  int colid=threadIdx.y + blockIdx.y * blockDim.y;
  if(rowid < nstruct && colid <natoms)
  {
    for(int k=0;k<xsplit;k++)
    {
      for(int l=0;l<ysplit;l++)
      {
        for(int m=0;m<zsplit;m++)
        {
        if(A[colid*3+rowid*natoms*3] >= xtick[k*2] && A[colid*3+rowid*natoms*3] < xtick[k*2+1] && A[colid*3+(rowid*natoms*3)+1] >= ytick[l*2] && A[colid*3+(rowid*natoms*3)+1] < ytick[l*2+1] && A[colid*3+(rowid*natoms*3)+2] >= ztick[m*2] && A[colid*3+(rowid*natoms*3)+2] < ztick[m*2+1] && exch[colid+rowid*natoms]==1)
          {
            atomicAdd(&density[m+l*zsplit+k*ysplit*zsplit],1);
	    break;
          }
           
        }
         
      }
      
    }
    }
}