//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"

__global__ void velocity_calc(float *A,float *B,float *C,int nstruct,int natoms,int xsplit,int ysplit, int zsplit, float *xtick,float *ytick, float *ztick)
{
  int rowid=threadIdx.x + blockIdx.x * blockDim.x;
  int colid=threadIdx.y + blockIdx.y * blockDim.y;
  float chkx,chky,chkz,velx,vely,velz,velmag;
  if(rowid < nstruct && colid <natoms)
  {
  chkx=A[colid*3+rowid*natoms*3];
  chky=A[colid*3+(rowid*natoms*3)+1];
  chkz=A[colid*3+(rowid*natoms*3)+2];
  velx=B[colid*3+rowid*natoms*3];
  vely=B[colid*3+(rowid*natoms*3)+1];
  velz=B[colid*3+(rowid*natoms*3)+2];
  velmag=sqrt(velx*velx+vely*vely+velz*velz);  
    for(int k=0;k<xsplit;k++)
    {
      for(int l=0;l<ysplit;l++)
      {
        for(int m=0;m<zsplit;m++)
        {
//if(A[colid*3+rowid*natoms*3] > xtick[k*2] && A[colid*3+rowid*natoms*3] < xtick[k*2+1] && A[colid*3+(rowid*natoms*3)+1] > ytick[l*2] && A[colid*3+(rowid*natoms*3)+1] < ytick[l*2+1] && A[colid*3+(rowid*natoms*3)+2] > ztick[m*2] && A[colid*3+(rowid*natoms*3)+2] < ztick[m*2+1])        
if(chkx > xtick[k*2] && chkx < xtick[k*2+1] && chky > ytick[l*2] && chky < ytick[l*2+1] && chkz > ztick[m*2] && chkz < ztick[m*2+1])
          {
            atomicAdd(&C[m+l*zsplit+k*ysplit*zsplit],velmag);
            break;
          }
        }
      }
    }
  }
}

