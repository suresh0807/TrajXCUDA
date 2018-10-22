//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




//filter for densities. buggy - donot use

#include "cudatools.cuh"

__global__ void density_filter(int *density, int *densityf, int mainsplit,int sub1split, int sub2split)
{
  int rowid=threadIdx.x + blockIdx.x * blockDim.x;
  int colid=threadIdx.y + blockIdx.y * blockDim.y;
  int widid=threadIdx.z + blockIdx.z * blockDim.z;
  int i1,i2,i3,j1,j2,j3,k1,k2,k3;
  if(rowid < mainsplit && widid < sub1split-1 && colid < sub2split-1 && colid != 0 && widid !=0)
  {
    i1=colid+widid*sub1split+rowid*sub2split*sub1split;
    j1=colid+(widid+1)*sub1split+rowid*sub2split*sub1split+1;
    k1=colid+(widid-1)*sub1split+rowid*sub2split*sub1split-1;
    i2=(colid-1)+widid*sub1split+rowid*sub2split*sub1split;
    j2=(colid-1)+(widid+1)*sub1split+rowid*sub2split*sub1split;
    k2=(colid-1)+(widid-1)*sub1split+rowid*sub2split*sub1split;
    i3=(colid+1)+widid*sub1split+rowid*sub2split*sub1split;
    j3=(colid+1)+(widid+1)*sub1split+rowid*sub2split*sub1split;
    k3=(colid+1)+(widid-1)*sub1split+rowid*sub2split*sub1split;
    densityf[colid+widid*mainsplit+rowid*sub2split*sub1split] = (density[i1]+density[i2]+density[i3]+density[j1]+density[j2]+density[j3]+density[k1]+density[k2]+density[k3])/9.0;
  }
}

