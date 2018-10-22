//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"

__global__ void SD_calc(float *A,float *B,float *C,int SD_store,int nElt1, int i, int j,int origins,int skips, int xsrt,int xend, int xski)
{
  int rowid=threadIdx.x + blockIdx.x * blockDim.x;
  //int colid=threadIdx.y + blockIdx.y * blockDim.y;
  int widid=threadIdx.z + blockIdx.z * blockDim.z;
  float chk,z;
  if(rowid < nElt1  &&  widid < origins/skips)
  {

    z=0.0;
    for(int k=xsrt; k< xend ; k=k+xski)
    //for(int k=0; k< 3 ; k=k+1)
    {
      chk=fabs(A[rowid*3+(widid*skips*nElt1*3)+(i*nElt1*3*SD_store)+k]-B[rowid*3+(j*nElt1*3)+(widid*skips*nElt1*3)+(i*nElt1*3*SD_store)+k]);
      z+=(chk * chk);
    }
//     printf("%f \n",C[rowid+(j*nElt1)+(widid*nElt1*SD_store)]);

    C[rowid+(j*nElt1)+widid*nElt1*SD_store] += z ;


  }
  
  
 /* int rowid=threadIdx.x + blockIdx.x * blockDim.x;
  int colid=threadIdx.y + blockIdx.y * blockDim.y;
  int widid=threadIdx.z + blockIdx.z * blockDim.z;
  float chk,z;
  if(rowid < nElt1  &&  rowid <SD_store && widid < origins/skips)
  {

    z=0.0;
    for(int k=0; k< 3 ; k++)
    {
      chk=fabs(A[rowid*3+(widid*skips*nElt1*3)+(i*nElt1*3*SD_store)+k]-B[rowid*3+(colid*nElt1*3)+(widid*skips*nElt1*3)+(i*nElt1*3*SD_store)+k]);
      z+=(chk * chk);
    }
//     printf("%f \n",C[rowid+(j*nElt1)+(widid*nElt1*SD_store)]);

    C[rowid+(colid*nElt1)+widid*nElt1*SD_store] += z ;


  }*/
  
}


__global__ void SD_calc(float *A,float *B,float *C,int *D, int *E, int *F, int SD_store,int nElt1, int i, int j,int origins,int skips, int xsrt,int xend, int xski,int whichwater)
{
  int rowid=threadIdx.x + blockIdx.x * blockDim.x;
  //int colid=threadIdx.y + blockIdx.y * blockDim.y;
  int widid=threadIdx.z + blockIdx.z * blockDim.z;
  float chk,z;
  if(rowid < nElt1  &&  widid < origins/skips)
  {
    if(D[rowid+(widid*skips*nElt1)+(i*nElt1*SD_store)] ==whichwater && D[rowid+(j*nElt1)+(widid*skips*nElt1)+(i*nElt1*SD_store)] ==whichwater && F[rowid+widid*nElt1] == 1)
    {
    z=0.0;
    for(int k=xsrt; k< xend ; k=k+xski)
    {
      chk=fabs(A[rowid*3+(widid*skips*nElt1*3)+(i*nElt1*3*SD_store)+k]-B[rowid*3+(j*nElt1*3)+(widid*skips*nElt1*3)+(i*nElt1*3*SD_store)+k]);
      z+=(chk * chk);
    }
    C[rowid+(j*nElt1)+widid*nElt1*SD_store] += z ;
    atomicAdd(&E[j+widid*SD_store],1);
    } 
    
    else if(whichwater==2)
    {
    z=0.0;
    for(int k=xsrt; k< xend ; k=k+xski)
    {
      chk=fabs(A[rowid*3+(widid*skips*nElt1*3)+(i*nElt1*3*SD_store)+k]-B[rowid*3+(j*nElt1*3)+(widid*skips*nElt1*3)+(i*nElt1*3*SD_store)+k]);
      z+=(chk * chk);
    }
     C[rowid+(j*nElt1)+widid*nElt1*SD_store] += z ; 
     atomicAdd(&E[j+widid*SD_store],1); // counting nelt1 in different time origins
    }
    
    else
    {
     F[rowid+widid*nElt1] = 0; 
    }
  }
}

 __global__ void SD_calc(float *A,float *B,float *C,int *D, int *E, int *F, int SD_store,int nElt1, int i, int j,int origins,int skips, int xsrt,int xend, int xski,int whichwater,int chase)
{
  int rowid=threadIdx.x + blockIdx.x * blockDim.x;
  //int colid=threadIdx.y + blockIdx.y * blockDim.y;
  int widid=threadIdx.z + blockIdx.z * blockDim.z;
  float chk,z;
  if(rowid < nElt1  &&  widid < origins/skips)
  {
    if(D[rowid+(widid*skips*nElt1)+(i*nElt1*SD_store)] ==whichwater && F[rowid+widid*nElt1] == 1)
    {
    z=0.0;
    for(int k=xsrt; k< xend ; k=k+xski)
    {
      chk=fabs(A[rowid*3+(widid*skips*nElt1*3)+(i*nElt1*3*SD_store)+k]-B[rowid*3+(j*nElt1*3)+(widid*skips*nElt1*3)+(i*nElt1*3*SD_store)+k]);
      z+=(chk * chk);
    }
    C[rowid+(j*nElt1)+widid*nElt1*SD_store] += z ;
    atomicAdd(&E[j+widid*SD_store],1);
    } 
    
    else if(whichwater==2)
    {
    z=0.0;
    for(int k=xsrt; k< xend ; k=k+xski)
    {
      chk=fabs(A[rowid*3+(widid*skips*nElt1*3)+(i*nElt1*3*SD_store)+k]-B[rowid*3+(j*nElt1*3)+(widid*skips*nElt1*3)+(i*nElt1*3*SD_store)+k]);
      z+=(chk * chk);
    }
     C[rowid+(j*nElt1)+widid*nElt1*SD_store] += z ; 
     atomicAdd(&E[j+widid*SD_store],1); // counting nelt1 in different time origins
    }
    
    else
    {
     F[rowid+widid*nElt1] = 0; 
    }
  } 
 /* int rowid=threadIdx.x + blockIdx.x * blockDim.x;
  int colid=threadIdx.y + blockIdx.y * blockDim.y;
  int widid=threadIdx.z + blockIdx.z * blockDim.z;
  float chk,z;
  if(rowid < nElt1  &&  rowid <SD_store && widid < origins/skips)
  {

    z=0.0;
    for(int k=0; k< 3 ; k++)
    {
      chk=fabs(A[rowid*3+(widid*skips*nElt1*3)+(i*nElt1*3*SD_store)+k]-B[rowid*3+(colid*nElt1*3)+(widid*skips*nElt1*3)+(i*nElt1*3*SD_store)+k]);
      z+=(chk * chk);
    }
//     printf("%f \n",C[rowid+(j*nElt1)+(widid*nElt1*SD_store)]);

    C[rowid+(colid*nElt1)+widid*nElt1*SD_store] += z ;


  }*/
  
}










__global__ void SD_calc(float *A,float *B,float *C,int SD_store,int nElt1, int i,int nstruct,int skips, int xsrt,int xend, int xski)
{
  int rowid=threadIdx.x + blockIdx.x * blockDim.x;
  int widid=threadIdx.z + blockIdx.z * blockDim.z;
  float chk,z;
  if(rowid < nElt1  &&  widid < nstruct)
  {
    z=0.0;
    for(int k=xsrt; k< xend ; k=k+xski)
    //for(int k=0; k< 3 ; k=k+1)
    {
      chk=fabs(A[rowid*3+(widid*nElt1*3)+k]-B[rowid*3+(widid*nElt1*3)+(i*nElt1*3)+k]);
      z+=(chk * chk);
    }
//     printf("%f \n",C[rowid+(j*nElt1)+(widid*nElt1*SD_store)]);

    C[rowid+widid*nElt1] = z ;


  }
}


