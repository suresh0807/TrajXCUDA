//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"

__global__ void Hbondmat(float *x, float*y, int *a, int *b, float bondist, int nElt1, int nElt2, int nstruct, float xvec, float yvec, float zvec)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int colid=threadIdx.y + blockIdx.y*blockDim.y;
  int widid=threadIdx.z + blockIdx.z*blockDim.z;
  //float vec[3],vecall;
  float chk;
  float z;
  //int chker;
  if(rowid < nElt1 && colid < nElt2 && widid < nstruct)
  {
    z=0.0;
    for(int k=0; k< 3 ; k++)
    {
      chk=fabs(x[rowid*3+widid*nElt1*3+k]-y[colid*3+widid*nElt2*3+k]);
      //vecall=x[rowid*3+widid*nElt1*3+k]-y[colid*3+widid*nElt2*3+k];
      if(k==0) {if(chk > xvec/2) {chk = chk - xvec; 
	                        //  if(vecall >0){vec[0]=y[colid*3+widid*nElt2*3+k]+ xvec;}
	                        //  else {vec[0]=y[colid*3+widid*nElt2*3+k]- xvec;}}}
      }}
      if(k==1) {if(chk > yvec/2) {chk = chk - yvec; 
	//                          if(vecall >0){vec[1]=y[colid*3+widid*nElt2*3+k]+ yvec;}
	//                          else {vec[1]=y[colid*3+widid*nElt2*3+k]- yvec;}}}
      }}
      if(k==2) {if(chk > zvec/2) {chk = chk - zvec; 
	//                          if(vecall >0){vec[2]=y[colid*3+widid*nElt2*3+k]+ zvec;}
	//                          else {vec[2]=y[colid*3+widid*nElt2*3+k]- zvec;}}}
      }}
      z+=(chk * chk);
    }
    z=sqrt(z);
   // printf(" %d %d %d %f %f \n",rowid,colid,widid,z,bondist);
    if (z < bondist)
    {
       //chker = a[rowid+widid*nElt1];
        //if(a[rowid+widid*nElt1] == 0) {b[rowid*2+widid*nElt1*2]=colid;b[rowid*2+widid*nElt1*2+1]=colid;}
        
       //b[rowid*2+widid*nElt1*2+1]=colid;
       //if(widid ==4998){chker = a[rowid+widid*nElt1]; printf("%d %d %d %d %f\n",widid,rowid,chker,colid,z);}
       //a[rowid+widid*nElt1]++;
       atomicAdd(&a[rowid+widid*nElt1],1);
       atomicCAS(&b[rowid*2+widid*nElt1*2],-1,colid);
       atomicCAS(&b[rowid*2+widid*nElt1*2+1],-1,colid);
       atomicCAS(&b[rowid*2+widid*nElt1*2+1],b[rowid*2+widid*nElt1*2],colid);
  
       //if(a[rowid+widid*nElt1] != 0) {b[rowid*2+widid*nElt1*2+1]=colid;}
       //if(widid ==4998){printf("%d %d %f\n",rowid,colid,z);}
           //for(int k=0; k< 2 ; k++){
           //{b[rowid*6+widid*nElt1*6+a[rowid+widid*nElt1]*3+k]=y[colid*3+widid*nElt2*3+k];}
           //{ if(a[rowid+widid*nElt1] == 1){b[rowid*6+widid*nElt1*6+k]=vec[k];}
	   //  else if(a[rowid+widid*nElt1] == 2){b[rowid*6+widid*nElt1*6+3+k]=vec[k];} }

	   //                        }
     }
    //     if (z < bondist)
    //{
    // if(widid ==4998){chker = a[rowid+widid*nElt1]; printf("%d %d %d %d %f\n",widid,rowid,chker,colid,z);}
      //chker = a[rowid+widid*nElt1];
      //b[rowid*2+widid*nElt1*2+chker-1]=colid;
    //}
      
    }
}

__global__ void Hbondmat(float *x, float*y, int *a, int *b, float *c,float bondist, int nElt1, int nElt2, int nstruct, float *latti)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int colid=threadIdx.y + blockIdx.y*blockDim.y;
  int widid=threadIdx.z + blockIdx.z*blockDim.z;
  float chk,xvec,yvec,zvec;
  float z;
  if(rowid < nElt1 && colid < nElt2 && widid < nstruct)
  {
    z=0.0;
    chk=0.0;
    xvec=latti[widid*6];
    yvec=latti[1+widid*6];
    zvec=latti[2+widid*6];
    for(int k=0; k< 3 ; k++)
    {
      chk=fabs(x[rowid*3+widid*nElt1*3+k]-y[colid*3+widid*nElt2*3+k]);
      if(k==0) {if(chk > xvec/2) {chk = chk - xvec;}}
      else if(k==1) {if(chk > yvec/2) {chk = chk - yvec;}}
      else if(k==2) {if(chk > zvec/2) {chk = chk - zvec;}}
      z+=(chk * chk);
    }
    z=sqrt(z);
    
    // I cannot save the distances directly because, 
    // I am unable to control the flow, I see a zero sometimes :(
    if (z < bondist)
    {
       atomicAdd(&a[rowid+widid*nElt1],1);
       /*
       atomicCAS(&b[rowid*4+widid*nElt1*4],-1,colid);
       atomicCAS(&b[rowid*4+widid*nElt1*4+1],-1,colid);
       atomicCAS(&b[rowid*4+widid*nElt1*4+2],-1,colid);
       atomicCAS(&b[rowid*4+widid*nElt1*4+3],-1,colid);
       else if(a[rowid+widid*nElt1]==2)
       {
       atomicCAS(&b[rowid*4+widid*nElt1*4+1],b[rowid*4+widid*nElt1*4],colid);
       }
       else if(a[rowid+widid*nElt1]==3)
       {
       atomicCAS(&b[rowid*4+widid*nElt1*4+2],b[rowid*4+widid*nElt1*4],colid);
       }
       else if(a[rowid+widid*nElt1]==4)
       {
       atomicCAS(&b[rowid*4+widid*nElt1*4+3],b[rowid*4+widid*nElt1*4],colid);
       }
       */
     }
    }
}

__global__ void Hbondmat(float *x, float*y, float *c,float bondist, int nElt1, int nElt2, float *latti,int widid)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int colid=threadIdx.y + blockIdx.y*blockDim.y;
  //int widid=threadIdx.z + blockIdx.z*blockDim.z;
  float chk,xvec,yvec,zvec;
  float z;
  if(rowid < nElt1 && colid < nElt2 )
  {
    z=0.0;
    chk=0.0;
    xvec=latti[widid*6];
    yvec=latti[1+widid*6];
    zvec=latti[2+widid*6];
    for(int k=0; k< 3 ; k++)
    {
      chk=fabs(x[rowid*3+widid*nElt1*3+k]-y[colid*3+widid*nElt2*3+k]);
      if(k==0)      {if(chk > xvec/2) {chk = chk - xvec;}}
      else if(k==1) {if(chk > yvec/2) {chk = chk - yvec;}}
      else if(k==2) {if(chk > zvec/2) {chk = chk - zvec;}}
      z+=(chk * chk);
    }
    z=sqrt(z);
    
    if (z < bondist)
    {
       c[colid+rowid*nElt2] = z;
     }
    }
}

//Atop orientation
__global__ void Hbondmat(float *x, float*y, float *c,float bondist, int nElt1, int nElt2, float *latti,int widid,int dirn)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int colid=threadIdx.y + blockIdx.y*blockDim.y;
  //int widid=threadIdx.z + blockIdx.z*blockDim.z;
  float chk,xvec,yvec,zvec,chkz;
  float z;
  if(rowid < nElt1 && colid < nElt2 )
  {
    z=0.0;
    chk=0.0;
    chkz=0.0;
    xvec=latti[widid*6];
    yvec=latti[1+widid*6];
    zvec=latti[2+widid*6];
    chkz=fabs(x[rowid*3+widid*nElt1*3+dirn]-y[colid*3+widid*nElt2*3+dirn]);
    if (chkz <= bondist)
    {
    for(int k=0; k< 3 ; k++)
    {
      chk=fabs(x[rowid*3+widid*nElt1*3+k]-y[colid*3+widid*nElt2*3+k]);
      if(k==0)      {if(chk > xvec/2) {chk = chk - xvec;}}
      else if(k==1) {if(chk > yvec/2) {chk = chk - yvec;}}
      else if(k==2) {if(chk > zvec/2) {chk = chk - zvec;}}
      z+=(chk * chk);
    }
    z=sqrt(z);
    c[colid+rowid*nElt2] = z;
    
     }
    }
}


__global__ void Hbondmat(float *x, float*y, int*exch, float *c,float bondist, int nElt1, int nElt2, float *latti,int widid)
{
  int rowid=threadIdx.x + blockIdx.x*blockDim.x;
  int colid=threadIdx.y + blockIdx.y*blockDim.y;
  //int widid=threadIdx.z + blockIdx.z*blockDim.z;
  float chk,xvec,yvec,zvec;
  float z;
  if(rowid < nElt1 && colid < nElt2 &&exch[rowid+widid*nElt1] >0 )
  {
    z=0.0;
    chk=0.0;
    xvec=latti[widid*6];
    yvec=latti[1+widid*6];
    zvec=latti[2+widid*6];
    for(int k=0; k< 3 ; k++)
    {
      chk=fabs(x[rowid*3+widid*nElt1*3+k]-y[colid*3+widid*nElt2*3+k]);
      if(k==0) {if(chk > xvec/2) {chk = chk - xvec;}}
      else if(k==1) {if(chk > yvec/2) {chk = chk - yvec;}}
      else if(k==2) {if(chk > zvec/2) {chk = chk - zvec;}}
      z+=(chk * chk);
    }
    z=sqrt(z);
    
    if (z < bondist)
    {
       c[colid+rowid*nElt2] = z;
     }
    }
}