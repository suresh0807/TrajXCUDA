//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"

__global__ void covermat(float *x, float*y, int *a, float bondist, int nElt1, int nElt2, int nstruct, float xvec, float yvec, float zvec)
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
    //if(rowid ==0){
    //printf(" %d %d %d %f %f \n",rowid,colid,widid,z,bondist);}
    if (z < bondist+0.01 )
    {
      atomicAdd(&a[rowid+widid*nElt1],1);
    }

  /*if(y[colid*3+widid*nElt2*3] > x[rowid*3+widid*nElt1*3] && y[colid*3+widid*nElt2*3+1] > x[rowid*3+widid*nElt1*3+1] && y[colid*3+widid*nElt2*3+2] > x[rowid*3+widid*nElt1*3+2] )
  {
   atomicAdd(&a[rowid+widid*nElt1],1);
  }
  else if(y[colid*3+widid*nElt2*3] < x[rowid*3+widid*nElt1*3] && y[colid*3+widid*nElt2*3+1] < x[rowid*3+widid*nElt1*3+1] && y[colid*3+widid*nElt2*3+2] < x[rowid*3+widid*nElt1*3+2] )
  {
   atomicAdd(&a[rowid+widid*nElt1],1);
  }
 */
 }
}

__global__ void covermat(float *x, float*y, int *a, float bondist, int nElt1, int nElt2, int nstruct, float *latti)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int colid=threadIdx.y + blockIdx.y*blockDim.y;
  int widid=threadIdx.z + blockIdx.z*blockDim.z;
  float chk,xvec,yvec,zvec;
  float z;
  if(rowid < nElt1 && colid < nElt2 && widid < nstruct)
  {
    z=0.0;
    xvec=latti[widid*6];
    yvec=latti[1+widid*6];
    zvec=latti[2+widid*6];
    for(int k=0; k< 3 ; k++)
    {
      chk=fabs(x[rowid*3+widid*nElt1*3+k]-y[colid*3+widid*nElt2*3+k]);
      if(k==0) {if(chk > xvec/2) {chk = chk - xvec;}}
      if(k==1) {if(chk > yvec/2) {chk = chk - yvec;}}
      if(k==2) {if(chk > zvec/2) {chk = chk - zvec;}}
      z+=(chk * chk);
    }
    z=sqrt(z);
    //if(rowid ==0){
    //printf(" %d %d %d %f %f \n",rowid,colid,widid,z,bondist);}
    if (z < bondist )
    {
      atomicAdd(&a[rowid+widid*nElt1],1);
    }

  /*if(y[colid*3+widid*nElt2*3] > x[rowid*3+widid*nElt1*3] && y[colid*3+widid*nElt2*3+1] > x[rowid*3+widid*nElt1*3+1] && y[colid*3+widid*nElt2*3+2] > x[rowid*3+widid*nElt1*3+2] )
  {
   atomicAdd(&a[rowid+widid*nElt1],1);
  }
  else if(y[colid*3+widid*nElt2*3] < x[rowid*3+widid*nElt1*3] && y[colid*3+widid*nElt2*3+1] < x[rowid*3+widid*nElt1*3+1] && y[colid*3+widid*nElt2*3+2] < x[rowid*3+widid*nElt1*3+2] )
  {
   atomicAdd(&a[rowid+widid*nElt1],1);
  }
 */
 }
}


__global__ void covermat(float *x, float*y, int *a, float minbondist, float maxbondist, int nElt1, int nElt2, int nstruct, float *latti, int dirn, float surf,float mintop, float maxtop)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int colid=threadIdx.y + blockIdx.y*blockDim.y;
  int widid=threadIdx.z + blockIdx.z*blockDim.z;
  float chk,xvec,yvec,zvec;
  float z;
  if(rowid < nElt1 && colid < nElt2 && widid < nstruct && x[rowid*3+widid*nElt1*3+dirn] >=mintop && x[rowid*3+widid*nElt1*3+dirn] <=maxtop&& y[colid*3+widid*nElt2*3+dirn] >=mintop && y[colid*3+widid*nElt2*3+dirn] <=maxtop)
  {
    z=0.0;
    xvec=latti[widid*6];
    yvec=latti[1+widid*6];
    zvec=latti[2+widid*6];
    for(int k=0; k< 3 ; k++)
    {
      chk=fabs(x[rowid*3+widid*nElt1*3+k]-y[colid*3+widid*nElt2*3+k]);
      if(k==0) {if(chk > xvec/2) {chk = chk - xvec;}}
      if(k==1) {if(chk > yvec/2) {chk = chk - yvec;}}
      if(k==2) {if(chk > zvec/2) {chk = chk - zvec;}}
      z+=(chk * chk);
    }
    z=sqrt(z);
    if ( z >= minbondist && z <= maxbondist &&x[rowid*3+widid*nElt1*3+dirn] > surf+minbondist && x[rowid*3+widid*nElt1*3+dirn] < surf+maxbondist)
    {
       atomicAdd(&a[rowid+widid*nElt1],1);
    }
 }
}
__global__ void covermat(float *x, float*y, int *a, float minbondist, float maxbondist, int nElt1, int nElt2, int nstruct, float *latti, float mintop, float maxtop,int dirn, float surf)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int colid=threadIdx.y + blockIdx.y*blockDim.y;
  int widid=threadIdx.z + blockIdx.z*blockDim.z;
  float chk,xvec,yvec,zvec;
  float z;
  if(rowid < nElt1 && colid < nElt2 && widid < nstruct&& y[colid*3+widid*nElt2*3+dirn] >mintop && y[colid*3+widid*nElt2*3+dirn] <maxtop && x[rowid*3+widid*nElt1*3+dirn] > maxtop)
  {
    z=0.0;
    xvec=latti[widid*6];
    yvec=latti[1+widid*6];
    zvec=latti[2+widid*6];
    for(int k=0; k< 3 ; k++)
    {
      chk=fabs(x[rowid*3+widid*nElt1*3+k]-y[colid*3+widid*nElt2*3+k]);
      if(k==0) {if(chk > xvec/2) {chk = chk - xvec;}}
      if(k==1) {if(chk > yvec/2) {chk = chk - yvec;}}
      if(k==2) {if(chk > zvec/2) {chk = chk - zvec;}}
      z+=(chk * chk);
    }
    z=sqrt(z);
    if (z > minbondist && z < maxbondist && x[rowid*3+widid*nElt1*3+dirn] > surf+minbondist && x[rowid*3+widid*nElt1*3+dirn] < surf+maxbondist)
    {
       atomicAdd(&a[rowid+widid*nElt1],1);
    }
 }
}


__global__ void covermat(float *x, int *a, float minbondist, float maxbondist, int nElt1, int nstruct, int dirn, float surf,float mintop, float maxtop)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int widid=threadIdx.z + blockIdx.z*blockDim.z;
  if(rowid < nElt1  && widid < nstruct && x[rowid*3+widid*nElt1*3+dirn] >mintop && x[rowid*3+widid*nElt1*3+dirn] <maxtop)
  {
    if (x[rowid*3+widid*nElt1*3+dirn] > surf+minbondist && x[rowid*3+widid*nElt1*3+dirn] < surf+maxbondist)
    {
       atomicAdd(&a[rowid+widid*nElt1],1);
    }
 }
}



__global__ void covermat(float *x, int *a, int nElt1, int nstruct, float minaz, float maxaz, int dirn)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int widid=threadIdx.z + blockIdx.z*blockDim.z;
  if(rowid < nElt1 && widid < nstruct)
  {
    if (x[rowid*3+widid*nElt1*3+dirn] < maxaz &&  x[rowid*3+widid*nElt1*3+dirn] > minaz)
    {
      atomicAdd(&a[rowid+widid*nElt1],1);
    }
 }
}

//residenttime
/*
__global__ void covermat(float *x, float*y, int *a, float bondist, int nElt1, int nElt2, int nstruct,float *latti,int dirn)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int colid=threadIdx.y + blockIdx.y*blockDim.y;
  int widid=threadIdx.z + blockIdx.z*blockDim.z;
  float chk,xvec,yvec,zvec,chkz;
  float z;
  if(rowid < nElt1 && colid < nElt2 && widid < nstruct)
  {
    z=0.0;
    xvec=latti[widid*6];
    yvec=latti[1+widid*6];
    zvec=latti[2+widid*6];
    chkz=fabs(x[rowid*3+widid*nElt1*3+dirn]-y[colid*3+widid*nElt2*3+dirn]);
    if (chkz <= bondist)
    {
    for(int k=0; k< 3 ; k++)
    {
      chk=fabs(x[rowid*3+widid*nElt1*3+k]-y[widid*3+k]);
      if(k==0) {if(chk > xvec/2) {chk = chk - xvec;}}
      else if(k==1) {if(chk > yvec/2) {chk = chk - yvec;}}
      else if(k==2) {if(chk > zvec/2) {chk = chk - zvec;}}
      z+=(chk * chk);
    }
    z=sqrt(z);
    
    if (z <= bondist+1.0 ) //(1.0 is tentative to allow for atoms between the spheres)
    {
      atomicAdd(&a[rowid+widid*nElt1],1);
    }
    }
 }
}
*/
__global__ void covermat(float *x, float*y, int *a, float bondist, int nElt1, int nElt2, int nstruct,float *latti,int dirn)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int colid=threadIdx.y + blockIdx.y*blockDim.y;
  int widid=threadIdx.z + blockIdx.z*blockDim.z;
  float chkz;
  if(rowid < nElt1 && colid < nElt2 && widid < nstruct)
  {
    chkz=fabs(x[rowid*3+widid*nElt1*3+dirn]-y[colid*3+widid*nElt2*3+dirn]);
    if(chkz <= bondist)  { atomicAdd(&a[rowid+widid*nElt1],1);}
  }
}

__global__ void covermat(float *x, float *y, int *a, float minbondist, float maxbondist, int index,float *latti, int nElt1, int nstruct)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int widid=threadIdx.z + blockIdx.z*blockDim.z;
  float chk,xvec,yvec,zvec;
  float z;
  if(rowid < nElt1 && widid < nstruct)
  {
    z=0.0;
    xvec=latti[widid*6];
    yvec=latti[1+widid*6];
    zvec=latti[2+widid*6];
    for(int k=0; k< 3 ; k++)
    {
      chk=fabs(x[rowid*3+widid*nElt1*3+k]-y[widid*3+k]);
      if(k==0) {if(chk > xvec/2) {chk = chk - xvec;}}
      if(k==1) {if(chk > yvec/2) {chk = chk - yvec;}}
      if(k==2) {if(chk > zvec/2) {chk = chk - zvec;}}
      z+=(chk * chk);
    }
    z=sqrt(z);
    if (z >= minbondist  &&  z <= maxbondist )
    {
      atomicAdd(&a[rowid+widid*nElt1],1);
    }
  }
}
