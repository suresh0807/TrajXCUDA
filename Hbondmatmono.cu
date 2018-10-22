//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"

__global__ void Hbondmatmono(float *x, float*y, int *a, int *b, float bondist, int nElt1, int nElt2, int nstruct, float xvec, float yvec, float zvec, float yxvec)
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
    z=sqrt((chkx*chkx) + (chky*chky) + (chkz*chkz));
    if (z < bondist)
    {
       atomicAdd(&a[rowid+widid*nElt1],1);
       atomicCAS(&b[rowid*2+widid*nElt1*2],-1,colid);
       atomicCAS(&b[rowid*2+widid*nElt1*2+1],-1,colid);
       atomicCAS(&b[rowid*2+widid*nElt1*2+1],b[rowid*2+widid*nElt1*2],colid);
     }
    }
}

__global__ void Hbondmatmono(float *x, float*y, int *a, int *b, float bondist, int nElt1, int nElt2, int nstruct, float xvec, float yvec, float zvec, float yxvec, float *latti)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int colid=threadIdx.y + blockIdx.y*blockDim.y;
  int widid=threadIdx.z + blockIdx.z*blockDim.z;
  float chkx,chky,chkz;
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
    z=sqrt((chkx*chkx) + (chky*chky) + (chkz*chkz));
    if (z < bondist)
    {
       atomicAdd(&a[rowid+widid*nElt1],1);
       atomicCAS(&b[rowid*2+widid*nElt1*2],-1,colid);
       atomicCAS(&b[rowid*2+widid*nElt1*2+1],-1,colid);
       atomicCAS(&b[rowid*2+widid*nElt1*2+1],b[rowid*2+widid*nElt1*2],colid);
     }
    }
}

__global__ void Hbondmatmono(float *x, float*y, float *c,float bondist, int nElt1, int nElt2, float *latti,int widid)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int colid=threadIdx.y + blockIdx.y*blockDim.y;
  //int widid=threadIdx.z + blockIdx.z*blockDim.z;
  float xvec,yvec,zvec,yxvec;
  float chkx,chky,chkz;
  float z;
  if(rowid < nElt1 && colid < nElt2 )
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
    z=sqrt((chkx*chkx) + (chky*chky) + (chkz*chkz));
    if (z < bondist)
    {
       c[colid+rowid*nElt2] = z;
    }
   }
}

__global__ void Hbondmatmono(float *x, float*y, float *c,float bondist, int nElt1, int nElt2, float *latti,int widid,int dirn)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int colid=threadIdx.y + blockIdx.y*blockDim.y;
  //int widid=threadIdx.z + blockIdx.z*blockDim.z;
  float xvec,yvec,zvec,yxvec,chkrr;
  float chkx,chky,chkz;
  float z;
  if(rowid < nElt1 && colid < nElt2 )
  {
    z=0.0;chkrr=0.0;
    xvec=latti[widid*6];
    yvec=latti[1+widid*6];
    zvec=latti[2+widid*6];
    yxvec=latti[3+widid*6];
    chkrr=fabs(x[rowid*3+widid*nElt1*3+dirn]-y[colid*3+widid*nElt2*3+dirn]);
    if (chkrr <= bondist)
    {
    chkx=fabs(x[rowid*3+widid*nElt1*3]-y[colid*3+widid*nElt2*3]);
    chky=fabs(x[rowid*3+widid*nElt1*3+1]-y[colid*3+widid*nElt2*3+1]);
    chkz=fabs(x[rowid*3+widid*nElt1*3+2]-y[colid*3+widid*nElt2*3+2]);
    if(chky > yvec/2) {chky = chky - yvec; chkx = chkx - yxvec;}
    if(chkx > xvec/2) {chkx = chkx - xvec;}
    if(chkz > zvec/2) {chkz = chkz - zvec;}
    z=sqrt((chkx*chkx) + (chky*chky) + (chkz*chkz));
    
    c[colid+rowid*nElt2] = z;
        }
   }
}


__global__ void Hbondmatmono(float *x, float*y, int*exch, float *c,float bondist, int nElt1, int nElt2, float *latti,int widid)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int colid=threadIdx.y + blockIdx.y*blockDim.y;
  //int widid=threadIdx.z + blockIdx.z*blockDim.z;
  float xvec,yvec,zvec,yxvec;
  float chkx,chky,chkz;
  float z;
  if(rowid < nElt1 && colid < nElt2 &&exch[rowid+widid*nElt1]>0)
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
    z=sqrt((chkx*chkx) + (chky*chky) + (chkz*chkz));
   
    if (z < bondist)
    {
       c[colid+rowid*nElt2] = z;
     }
   }
}