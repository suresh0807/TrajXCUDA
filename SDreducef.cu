//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"

__global__ void SDreducef(float *a, float *b, int nElt1, int SD_store, int origins,int skips)
{
  int rowid=threadIdx.x + blockIdx.x * blockDim.x;
  int colid=threadIdx.y + blockIdx.y * blockDim.y;
  if(rowid < nElt1 && colid < SD_store)
  {
    for(int i = 0; i < origins/skips; i++)
    {
      //printf("%f \n", a[rowid+(colid*nElt1)+(i*nElt1*SD_store)]); 
      b[rowid+colid*nElt1] += a[rowid+(colid*nElt1)+(i*nElt1*SD_store)];
    }
    b[rowid+colid*nElt1] = b[rowid+colid*nElt1]/((float(origins)/float(skips)));
  }
}

__global__ void SDreducef(float *a, float *b, int nElt1, int SD_store, int nstruct,int skips,int i)
{
  int rowid=threadIdx.x + blockIdx.x * blockDim.x;
  if(rowid < nElt1 )
  {
    for(int i1 = 0; i1 < nstruct; i1++)
    {
      b[rowid+i*nElt1] += a[rowid+i1*nElt1];
    }
    b[rowid+i*nElt1] = b[rowid+i*nElt1]/float(nstruct);
  }
  


/* extern __shared__ int share[];
 unsigned int tid = threadIdx.x;
 unsigned int rowid = threadIdx.x + blockIdx.x * blockDim.x;
 share[tid] = a[rowid];
 __syncthreads();
 
 for(int i=1; i<blockDim.x ; i *=2)
 {
   if(tid % (2*i) ==0)
   {
     share[tid] += share[tid+i];
   }
   __syncthreads();
 }
 
 if(tid ==0) { b[blockIdx.x] = share[0];}
   if(rowid < nElt1 )
  {
   b[rowid+i*nElt1] = b[rowid+i*nElt1]/float(nstruct);
  }
*/

}

