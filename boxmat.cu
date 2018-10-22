//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"


__global__ void boxmat(float *x, int *a, int nElt1, int nstruct,float maxx, float minx, float maxy,float miny, float maxz, float minz)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int widid=threadIdx.z + blockIdx.z*blockDim.z;
  
  if(rowid < nElt1 && widid < nstruct)
  {

   if( x[rowid*3+widid*nElt1*3] < maxx && x[rowid*3+widid*nElt1*3] > minx &&  x[rowid*3+widid*nElt1*3+1] < maxy && x[rowid*3+widid*nElt1*3+1] > miny && x[rowid*3+widid*nElt1*3+2] < maxz && x[rowid*3+widid*nElt1*3+2] > minz )
   {
     atomicAdd(&a[rowid+widid*nElt1],1);
   }
   
   }
}