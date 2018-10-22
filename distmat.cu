//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"

__global__ void distmat(float *x, float*y, float *a, float bin_rad, int bin, int nElt1, int nElt2, int nstruct, float xvec, float yvec, float zvec, float yxvec, float rdf_max_rad)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int colid=threadIdx.y + blockIdx.y*blockDim.y;
  int widid=threadIdx.z + blockIdx.z*blockDim.z;
  float chk;
  float z;
  float maxval;
  int index;
  if(rowid < nElt1 && colid < nElt2 && widid < nstruct)
  {
    //maxval = xvec/2;
    maxval = rdf_max_rad;
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
    index = z * bin / maxval;
    for(int loop=index;loop<bin;loop++)
    {
      if(loop == index) {atomicAdd(&a[index+rowid*bin+widid*nElt1*bin],1.0);break;}
    }
  }
}


__global__ void distmat(float *x, float*y, float *a, int *exch1, int *exch2, float bin_rad, int bin, int nElt1, int nElt2, int nstruct, float xvec, float yvec, float zvec, float yxvec, float rdf_max_rad, int whichwater)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int colid=threadIdx.y + blockIdx.y*blockDim.y;
  int widid=threadIdx.z + blockIdx.z*blockDim.z;
  float chk;
  float z;
  float maxval;
  int index;
  if(rowid < nElt1 && colid < nElt2 && widid < nstruct)
  {
    //maxval = xvec/2;
    maxval = rdf_max_rad;
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
    index = z * bin / maxval;
    if(exch1[rowid+nElt1*widid]==whichwater && exch2[colid+nElt2*widid]==whichwater)
    {
    for(int loop=index;loop<bin;loop++)
    {
      if(loop == index) {atomicAdd(&a[index+rowid*bin+widid*nElt1*bin],1.0);break;}
    }
    }
    else if (whichwater ==2)
    {
    for(int loop=index;loop<bin;loop++)
    {
      if(loop == index) {atomicAdd(&a[index+rowid*bin+widid*nElt1*bin],1.0);break;}
    }
    }
  }
}

// Fork bulk region --- exch2 is not needed.
__global__ void distmat(float *x, float*y, float *a, int *exch1, int *exch2, float bin_rad, int bin, int nElt1, int nElt2, int nstruct, float *latti, float rdf_max_rad, int whichwater)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int colid=threadIdx.y + blockIdx.y*blockDim.y;
  int widid=threadIdx.z + blockIdx.z*blockDim.z;
  float chk,xvec,yvec,zvec;
  float z;
  float maxval;
  int index;
  if(rowid < nElt1 && colid < nElt2 && widid < nstruct)
  {
    //maxval = xvec/2;
    maxval = rdf_max_rad;
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
    index = z * bin / maxval;
    if(exch1[rowid+nElt1*widid]==whichwater)// && exch2[colid+nElt2*widid]==whichwater)
    {
    for(int loop=index;loop<bin;loop++)
    {
      if(loop == index) {atomicAdd(&a[index+rowid*bin+widid*nElt1*bin],1.0);break;}
    }
    }
    else if (whichwater ==2)
    {
    for(int loop=index;loop<bin;loop++)
    {
      if(loop == index) {atomicAdd(&a[index+rowid*bin+widid*nElt1*bin],1.0);break;}
    }
    }
  }
}

__global__ void distmat(float *x, float*y, float *a, int *exch1, int *exch2, float bin_rad, int bin, int nElt1, int nElt2, int nstruct, float *latti, float rdf_max_rad, int whichwater, float mid, int dirn)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int colid=threadIdx.y + blockIdx.y*blockDim.y;
  int widid=threadIdx.z + blockIdx.z*blockDim.z;
  float chk,xvec,yvec,zvec;
  float z;
  float maxval;
  int index;
  if(rowid < nElt1 && colid < nElt2 && widid < nstruct)
  {
    //maxval = xvec/2;
    maxval = rdf_max_rad;
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
    index = z * bin / maxval;
    if(exch1[rowid+nElt1*widid]==whichwater && x[rowid*3+widid*nElt1*3+dirn] <= mid  && y[colid*3+widid*nElt2*3+dirn] >= x[rowid*3+widid*nElt1*3+dirn] )//&& y[colid*3+widid*nElt2*3+2] <= mid)// && exch2[colid+nElt2*widid]==whichwater)
    {
    for(int loop=index;loop<bin;loop++)
    {
      if(loop == index) {atomicAdd(&a[index+rowid*bin+widid*nElt1*bin],1.0);break;}
    }
    }
    else if(exch1[rowid+nElt1*widid]==whichwater && x[rowid*3+widid*nElt1*3+dirn] >= mid  && y[colid*3+widid*nElt2*3+dirn] <= x[rowid*3+widid*nElt1*3+dirn])// && y[colid*3+widid*nElt2*3+2] >=mid && exch2[colid+nElt2*widid]==whichwater)
    {
    for(int loop=index;loop<bin;loop++)
    {
      if(loop == index) {atomicAdd(&a[index+rowid*bin+widid*nElt1*bin],1.0);break;}
    }
    }
    else if (whichwater ==2)
    {
    for(int loop=index;loop<bin;loop++)
    {
      if(loop == index) {atomicAdd(&a[index+rowid*bin+widid*nElt1*bin],1.0);break;}
    }
    }
    
  }
}

__global__ void distmat(float *x, float*y, float *a, int nElt1, int nElt2, int nstruct, float *latti, float cutoff)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int colid=threadIdx.y + blockIdx.y*blockDim.y;
  int widid=threadIdx.z + blockIdx.z*blockDim.z;
  float chk,xvec,yvec,zvec;
  float z;
  float index=0.0;
  if(rowid < nElt1 && colid < nElt2 && widid < nstruct && rowid != colid)
  {
    xvec=latti[widid*6];
    yvec=latti[1+widid*6];
    zvec=latti[2+widid*6];
    z=0.0;
    for(int k=0; k< 3 ; k++)
    {
      chk=fabs(x[rowid*4+widid*nElt1*4+k]-y[colid*4+widid*nElt2*4+k]);
      //if(rowid==1){printf("%f %f",x[rowid*4+widid*nElt1*4+k],y[colid*4+widid*nElt2*4+k]);}
      if(k==0) {if(chk > xvec/2) {chk = chk - xvec;}}
      if(k==1) {if(chk > yvec/2) {chk = chk - yvec;}}
      if(k==2) {if(chk > zvec/2) {chk = chk - zvec;}}
      //if(rowid==1){printf("%f %f",z,index);}
      z+=(chk * chk);
    }
    z=sqrt(z);
     if(z >0 && z < cutoff){index = (x[rowid*4+widid*nElt1*4+3]*y[colid*4+widid*nElt2*4+3]) / (z*1.88973);}
    if(rowid==1){printf("%f %f ",z,index);}
    a[colid+rowid*nElt2+widid*nElt1*nElt2]=index;
    
  }
}