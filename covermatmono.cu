//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"

__global__ void covermatmono(float *x, float*y, int *a, float bondist, int nElt1, int nElt2, int nstruct, float xvec, float yvec, float zvec, float yxvec)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int colid=threadIdx.y + blockIdx.y*blockDim.y;
  int widid=threadIdx.z + blockIdx.z*blockDim.z;
  float chkx,chky,chkz;
  float z;
  if(rowid < nElt1 && colid < nElt2 && widid < nstruct)
  {
    z=0.0;
    chkx=fabs(x[rowid*3+widid*nElt1*3]-y[colid*3+widid*nElt2*3]);
    chky=fabs(x[rowid*3+widid*nElt1*3+1]-y[colid*3+widid*nElt2*3+1]);
    chkz=fabs(x[rowid*3+widid*nElt1*3+2]-y[colid*3+widid*nElt2*3+2]);
    if(chky > yvec/2) {chky = chky - yvec; chkx = chkx - yxvec;}
    if(chkx > xvec/2) {chkx = chkx - xvec;}
    if(chkz > zvec/2) {chkz = chkz - zvec;}
    //z=sqrt(chk*chk);
    z=sqrt((chkx*chkx) + (chky*chky) + (chkz*chkz));
    //if(rowid ==2 && widid == 1){
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

__global__ void covermatmono(float *x, float*y, int *a, float bondist, int nElt1, int nElt2, int nstruct, float *latti)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int colid=threadIdx.y + blockIdx.y*blockDim.y;
  int widid=threadIdx.z + blockIdx.z*blockDim.z;
  float chkx,chky,chkz,xvec,yvec,zvec,yxvec;
  float z;
  if(rowid < nElt1 && colid < nElt2 && widid < nstruct)
  {
    z=0.0;
    xvec=latti[widid*6];
    yvec=latti[1+widid*6];
    zvec=latti[2+widid*6];
    yxvec=latti[3+widid*6];
    chkx=fabs(x[rowid*3+widid*nElt1*3]-y[colid*3+widid*nElt2*3]);
    chky=fabs(x[rowid*3+widid*nElt1*3+1]-y[colid*3+widid*nElt2*3+1]);
    chkz=fabs(x[rowid*3+widid*nElt1*3+2]-y[colid*3+widid*nElt2*3+2]);
    if(chky > yvec/2) {chky = chky - yvec; chkx = chkx - yxvec;}
    if(chkx > xvec/2) {chkx = chkx - xvec;}
    if(chkz > zvec/2) {chkz = chkz - zvec;}
    //z=sqrt(chk*chk);
    z=sqrt((chkx*chkx) + (chky*chky) + (chkz*chkz));
    //if(rowid ==2 && widid == 1){
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


__global__ void covermatmono(float *x, float*y, int *a, float minbondist, float maxbondist, int nElt1, int nElt2, int nstruct, float *latti , int dirn, float surf,float mintop, float maxtop  )
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int colid=threadIdx.y + blockIdx.y*blockDim.y;
  int widid=threadIdx.z + blockIdx.z*blockDim.z;
  float chkx,chky,chkz,xvec,yvec,zvec,yxvec;
  float z;
  if(rowid < nElt1 && colid < nElt2 && widid < nstruct && x[rowid*3+widid*nElt1*3+dirn] >mintop && x[rowid*3+widid*nElt1*3+dirn] <maxtop && y[colid*3+widid*nElt2*3+dirn] >mintop && y[colid*3+widid*nElt2*3+dirn] <maxtop)
  //if(rowid < nElt1 && colid < nElt2 && widid < nstruct)
  {
    z=0.0;
    xvec=latti[widid*6];
    yvec=latti[1+widid*6];
    zvec=latti[2+widid*6];
    yxvec=latti[3+widid*6];
    chkx=fabs(x[rowid*3+widid*nElt1*3]-y[colid*3+widid*nElt2*3]);
    chky=fabs(x[rowid*3+widid*nElt1*3+1]-y[colid*3+widid*nElt2*3+1]);
    chkz=fabs(x[rowid*3+widid*nElt1*3+2]-y[colid*3+widid*nElt2*3+2]);
    if(chky > yvec/2) {chky = chky - yvec; chkx = chkx - yxvec;}
    if(chkx > xvec/2) {chkx = chkx - xvec;}
    if(chkz > zvec/2) {chkz = chkz - zvec;}
    //z=sqrt(chk*chk);
    z=sqrt((chkx*chkx) + (chky*chky) + (chkz*chkz));
    //if(rowid ==2 && widid == 1){
    //printf(" %d %d %d %f %f \n",rowid,colid,widid,z,bondist);}
    if (z < maxbondist && z > minbondist && x[rowid*3+widid*nElt1*3+dirn] > surf+minbondist && x[rowid*3+widid*nElt1*3+dirn] < surf+maxbondist)
    //if (z < maxbondist && z > minbondist)
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

__global__ void covermatmono(float *x, float*y, int *a, float minbondist, float maxbondist, int nElt1, int nElt2, int nstruct, float *latti , float mintop, float maxtop,int dirn, float surf)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int colid=threadIdx.y + blockIdx.y*blockDim.y;
  int widid=threadIdx.z + blockIdx.z*blockDim.z;
  float chkx,chky,chkz,xvec,yvec,zvec,yxvec;
  float z;
  if(rowid < nElt1 && colid < nElt2 && widid < nstruct && widid < nstruct && y[colid*3+widid*nElt2*3+dirn] >mintop && y[colid*3+widid*nElt2*3+dirn] <maxtop&& x[rowid*3+widid*nElt1*3+dirn] > maxtop)
  //if(rowid < nElt1 && colid < nElt2 && widid < nstruct)
  {
    z=0.0;
    xvec=latti[widid*6];
    yvec=latti[1+widid*6];
    zvec=latti[2+widid*6];
    yxvec=latti[3+widid*6];
    chkx=fabs(x[rowid*3+widid*nElt1*3]-y[colid*3+widid*nElt2*3]);
    chky=fabs(x[rowid*3+widid*nElt1*3+1]-y[colid*3+widid*nElt2*3+1]);
    chkz=fabs(x[rowid*3+widid*nElt1*3+2]-y[colid*3+widid*nElt2*3+2]);
    if(chky > yvec/2) {chky = chky - yvec; chkx = chkx - yxvec;}
    if(chkx > xvec/2) {chkx = chkx - xvec;}
    if(chkz > zvec/2) {chkz = chkz - zvec;}
    //z=sqrt(chk*chk);
    z=sqrt((chkx*chkx) + (chky*chky) + (chkz*chkz));
    //if(rowid ==2 && widid == 1){
    //printf(" %d %d %d %f %f \n",rowid,colid,widid,z,bondist);}
    if (z < maxbondist && z > minbondist && x[rowid*3+widid*nElt1*3+dirn] > surf+minbondist && x[rowid*3+widid*nElt1*3+dirn] < surf+maxbondist)
    //if (z < maxbondist && z > minbondist)
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


__global__ void covermatmono(float *x, int *a, int nElt1, int nstruct, float minaz, float maxaz, int dirn)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int widid=threadIdx.z + blockIdx.z*blockDim.z;
  if(rowid < nElt1 &&  widid < nstruct)
  {
        if (x[rowid*3+widid*nElt1*3+dirn] < maxaz &&  x[rowid*3+widid*nElt1*3+dirn] > minaz)
    {
      atomicAdd(&a[rowid+widid*nElt1],1);
    }

 }
}

//residenttime
__global__ void covermatmono(float *x, float*y, int *a, float bondist, int nElt1, int nElt2, int nstruct,float *latti,int dirn)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int colid=threadIdx.y + blockIdx.y*blockDim.y;
  int widid=threadIdx.z + blockIdx.z*blockDim.z;
  float chkrr;
  
  if(rowid < nElt1 && colid < nElt2 && widid < nstruct)
  {
    chkrr=fabs(x[rowid*3+widid*nElt1*3+dirn]-y[colid*3+widid*nElt2*3+dirn]);
    if (chkrr <= bondist)
    {    
      atomicAdd(&a[rowid+widid*nElt1],1);
    }
 }
}

__global__ void covermatmono(float *x, float *y, int *a, float minbondist, float maxbondist, int index,float *latti, int nElt1, int nstruct)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int widid=threadIdx.z + blockIdx.z*blockDim.z;
  float chkx,chky,chkz,xvec,yvec,zvec,yxvec;
  float z;
  if(rowid < nElt1 && widid < nstruct)
  //if(rowid < nElt1 && colid < nElt2 && widid < nstruct)
  {
    z=0.0;
    xvec=latti[widid*6];
    yvec=latti[1+widid*6];
    zvec=latti[2+widid*6];
    yxvec=latti[3+widid*6];
    chkx=fabs(x[rowid*3+widid*nElt1*3]-y[widid*3]);
    chky=fabs(x[rowid*3+widid*nElt1*3+1]-y[widid*3+1]);
    chkz=fabs(x[rowid*3+widid*nElt1*3+2]-y[widid*3+2]);
    if(chky > yvec/2) {chky = chky - yvec; chkx = chkx - yxvec;}
    if(chkx > xvec/2) {chkx = chkx - xvec;}
    if(chkz > zvec/2) {chkz = chkz - zvec;}
    //z=sqrt(chk*chk);
    z=sqrt((chkx*chkx) + (chky*chky) + (chkz*chkz));
    //if(rowid ==2 && widid == 1){
    //printf(" %d %d %d %f %f \n",rowid,colid,widid,z,bondist);}
    if (z <= maxbondist && z >= minbondist)
    //if (z < maxbondist && z > minbondist)
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