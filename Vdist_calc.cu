//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################





#include "cudatools.cuh"

__global__ void Vdist_calc(float *A,float *B, float *density,float *xdensity,float *ydensity,float *zdensity,int nstruct,int natoms,int zsplit, float *ztick,int Dirn)
{
  int rowid=threadIdx.x + blockIdx.x * blockDim.x;
  int colid=threadIdx.y + blockIdx.y * blockDim.y;
  float vel;
  if(rowid < nstruct && colid <natoms)
  {
        for(int m=0;m<zsplit;m++)
        {
        if( A[colid*3+(rowid*natoms*3)+Dirn] >= ztick[m*2] && A[colid*3+(rowid*natoms*3)+Dirn] < ztick[m*2+1])
          {
	    vel= sqrt((B[colid*3+(rowid*natoms*3)]*B[colid*3+(rowid*natoms*3)]) + (B[colid*3+(rowid*natoms*3)+1]*B[colid*3+(rowid*natoms*3)+1]) + (B[colid*3+(rowid*natoms*3)+2]*B[colid*3+(rowid*natoms*3)+2]));
	    atomicAdd(&density[m],vel);
	    vel=fabs(B[colid*3+(rowid*natoms*3)]);
	    atomicAdd(&xdensity[m],vel);
	    vel=fabs(B[colid*3+(rowid*natoms*3)+1]);
	    atomicAdd(&ydensity[m],vel);
	    vel=fabs(B[colid*3+(rowid*natoms*3)+2]);
	    atomicAdd(&zdensity[m],vel);
            break;
          }
        }
   }
}