//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################





#include "cudatools.cuh"

__global__ void Hbond_calc(float *A,float *Hbond_density,float *angle_density,float *distance_density, float *OHdensity, float *O_Hdensity, float *angdensity, int nstruct,int natoms,int zsplit, float *ztick, float angle, float max_O_O,int Dirn)
{
  int rowid=threadIdx.x + blockIdx.x * blockDim.x;
  int colid=threadIdx.y + blockIdx.y * blockDim.y;
  if(rowid < nstruct && colid <natoms)
  {
        for(int m=0;m<zsplit;m++)
        {
        if( A[colid*9+(rowid*natoms*9)+Dirn] >= ztick[m*2] && A[colid*9+(rowid*natoms*9)+Dirn] < ztick[m*2+1])
          {
            atomicAdd(&Hbond_density[m],1.0);
            atomicAdd(&angle_density[m],A[colid*9+(rowid*natoms*9)+3]);
            atomicAdd(&distance_density[m],A[colid*9+(rowid*natoms*9)+4]);
	    atomicAdd(&OHdensity[m],A[colid*9+(rowid*natoms*9)+6]);
            atomicAdd(&O_Hdensity[m],A[colid*9+(rowid*natoms*9)+7]);
	    atomicAdd(&angdensity[m],A[colid*9+(rowid*natoms*9)+8]);
            break;
          }
        }
   }
}


