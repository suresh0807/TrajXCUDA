//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"

__global__ void Ofragmat(float *x, float*y, int *a, float bondist, int nElt1, int nElt2, int nstruct, float xvec, float yvec, float zvec)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int colid=threadIdx.y + blockIdx.y*blockDim.y;
  int widid=threadIdx.z + blockIdx.z*blockDim.z;
  float chk;
  float z;
  if(rowid < nElt1 && colid < nElt2 && widid < nstruct)
  {
    z=0.0;
    for(int k=0; k< 3 ; k++)
    {
      chk=fabs(x[rowid*3+widid*nElt1*3+k]-y[colid*3+widid*nElt2*3+k]);
      if(k==0) {if(chk > xvec/2) {chk = chk - xvec;}}
      if(k==1) {if(chk > yvec/2) {chk = chk - yvec;}}
      if(k==2) {if(chk > zvec/2) {chk = chk - zvec;}}
      z+=(chk * chk);
    }
    z=sqrt(z);
   // printf(" %d %d %d %f %f \n",rowid,colid,widid,z,bondist);
    if (z < bondist)
    {
    atomicAdd(&a[rowid+widid*nElt1],1);
    }
  }
}

