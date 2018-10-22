//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################



#include "cudatools.cuh"

void compute_densityvelocitygrid()
{
    if(plot=="int")        { nElt1 = nA_int;}
    else if(plot=="bulk")  { nElt1 = nA_bulk;}
    
cout << type <<" is chosen"<<endl;
cout<<"densityvelocitygrid: " << nElt1 << Elt1 << " atoms are there in each structure"<<endl;
    printf("%f %f %f %f %f %f\n",minx,maxx,miny,maxy,minz,maxz);
    
    if(set_max_z == "yes")
{
  minz = set_minz;
  maxz = set_maxz;
}
    if(set_max_x == "yes")
{
  minx = set_minx;
  maxx = set_maxx;
}
    if(set_max_y == "yes")
{
  miny = set_miny;
  maxy = set_maxy;
}
    xrange=maxx-minx;
    yrange=maxy-miny;
    zrange=maxz-minz;

    xint = xrange/xsplit;
    yint = yrange/ysplit;
    zint = zrange/zsplit;

    xtick=(float*) malloc (sizeof(float)*xsplit*2);
    ytick=(float*) malloc (sizeof(float)*ysplit*2);
    ztick=(float*) malloc (sizeof(float)*zsplit*2);
    for(int i=0;i<xsplit;i++)
    {
      xtick[i*2]=minx+(i*xint);
      xtick[i*2+1]=minx+((i+1)*xint);
    }
    for(int i=0;i<ysplit;i++)
    {
      ytick[i*2]=miny+(i*yint);
      ytick[i*2+1]=miny+((i+1)*yint);
    }
    for(int i=0;i<zsplit;i++)
    {
      ztick[i*2]=minz+(i*zint);
      ztick[i*2+1]=minz+((i+1)*zint);
    }

FILE *xyztick = fopen("xyztick.dat","wt");
    for(int i=0;i<xsplit;i++)
    {
      fprintf(xyztick,"%f %f \n",xtick[i*2],xtick[i*2+1]);
    }
      fprintf(xyztick," \n");
    for(int i=0;i<ysplit;i++)
    {
      fprintf(xyztick,"%f %f \n",ytick[i*2],ytick[i*2+1]);    
    }
      fprintf(xyztick," \n");
    for(int i=0;i<zsplit;i++)
    {
      fprintf(xyztick,"%f %f \n",ztick[i*2],ztick[i*2+1]);
    }
fclose(xyztick);

    density=(int *) malloc (sizeof(int)*xsplit*ysplit*zsplit);
    for(int i=0;i<xsplit;i++)
    {
      for(int j=0;j<ysplit;j++)
      {
        for(int k=0;k<zsplit;k++)
        {
          density[k+j*zsplit+i*ysplit*zsplit]=0;
        }
      }
    }
    
    velocity=(float *) malloc (sizeof(float)*xsplit*ysplit*zsplit);
    for(int i=0;i<xsplit;i++)
    {
      for(int j=0;j<ysplit;j++)
      {
        for(int k=0;k<zsplit;k++)
        {
          velocity[k+j*zsplit+i*ysplit*zsplit]=0.0;
        }
      }
    }
    
    
///*
    printf("Start of cuda calculation\n");
    cudaMalloc((void **)&dev_A,sizeof(float)*(nstruct/skip)*nElt1*3);
    cudaMalloc((void **)&dev_density,sizeof(int)*xsplit*ysplit*zsplit);
    cudaMalloc((void **)&dev_xtick, sizeof(float)*xsplit*2);
    cudaMalloc((void **)&dev_ytick, sizeof(float)*ysplit*2);
    cudaMalloc((void **)&dev_ztick, sizeof(float)*zsplit*2);
    if(plot=="int")
    {
      cudaMemcpy(dev_A,A_int,sizeof(float)*(nstruct/skip)*nA_int*3,cudaMemcpyHostToDevice);
    }
    else if(plot=="bulk")
    {
      cudaMemcpy(dev_A,A_bulk,sizeof(float)*(nstruct/skip)*nA_bulk*3,cudaMemcpyHostToDevice);
    }
    else if(plot=="all")
    {
      cudaMemcpy(dev_A,A,sizeof(float)*(nstruct/skip)*nElt1*3,cudaMemcpyHostToDevice);
    }
    cudaMemcpy(dev_density,density,sizeof(int)*xsplit*ysplit*zsplit,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_xtick,xtick,sizeof(float)*xsplit*2,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_ytick,ytick,sizeof(float)*ysplit*2,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_ztick,ztick,sizeof(float)*zsplit*2,cudaMemcpyHostToDevice);
    
    printf("%ld B needed\n",(((nstruct/skip)*nElt1*3+xsplit*ysplit*zsplit+xsplit*2+ysplit*2+zsplit*2)*4));
    dim3 dimBlock(32,32,1);
    //dim3 dimBlock(1,1,1);
    dim3 dimGrid(((nstruct/skip)+dimBlock.x-1)/dimBlock.x,(nElt1+dimBlock.y-1)/dimBlock.y,1);
    //dim3 dimGrid(nstruct,nElt1,1);
    density_calc<<<dimGrid,dimBlock>>>(dev_A,dev_density,(nstruct/skip),nElt1,xsplit,ysplit,zsplit,dev_xtick,dev_ytick,dev_ztick);
    cudaMemcpy(density,dev_density,sizeof(int)*xsplit*ysplit*zsplit,cudaMemcpyDeviceToHost);
    cudaFree(dev_density);
    
    
    
    cudaMalloc((void **)&dev_VEL,sizeof(float)*(nstruct/skip)*nElt1*3);
    cudaMalloc((void **)&dev_velocity,sizeof(float)*xsplit*ysplit*zsplit);
    cudaMemcpy(dev_VEL,VEL,sizeof(float)*(nstruct/skip)*nElt1*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_velocity,velocity,sizeof(float)*xsplit*ysplit*zsplit,cudaMemcpyHostToDevice);
    printf("%ld B needed\n",(((nstruct/skip)*nElt1*3+xsplit*ysplit*zsplit+xsplit*2+ysplit*2+zsplit*2)*4));
    velocity_calc<<<dimGrid,dimBlock>>>(dev_A,dev_VEL,dev_velocity,(nstruct/skip),nElt1,xsplit,ysplit,zsplit,dev_xtick,dev_ytick,dev_ztick);
    cudaMemcpy(velocity,dev_velocity,sizeof(float)*xsplit*ysplit*zsplit,cudaMemcpyDeviceToHost);
    cudaFree(dev_velocity);
    cudaFree(dev_VEL);
    cudaFree(dev_A);
    cudaFree(dev_xtick);
    cudaFree(dev_ytick);
    cudaFree(dev_ztick);
    
    
FILE *densityfile=fopen("out","wt");
FILE *velocityfile=fopen("velout","wt");

if(order=="XYZ") {

for(int l=0;l<ysplit;l++)
{
  //fprintf(densityfile,"%d %d\n",zsplit,xsplit);
  fprintf(densityfile,"#%d \n",l);
  fprintf(velocityfile,"#%d \n",l);
for(int m=0;m<zsplit;m++)
{
for(int k=0;k<xsplit;k++)
{
    fprintf(velocityfile,"%f ",velocity[m+l*zsplit+k*ysplit*zsplit]/(density[m+l*zsplit+k*ysplit*zsplit]+0.0001));
    fprintf(densityfile,"%d ",density[m+l*zsplit+k*ysplit*zsplit]);
}
fprintf(densityfile,"\n");
fprintf(velocityfile,"\n");
}
//fprintf(densityfile,"5\n");
//fprintf(densityfile,"1 10.\n");
//fprintf(densityfile,"2 15.\n");
//fprintf(densityfile,"3 20.\n");
//fprintf(densityfile,"4 25.\n");
//fprintf(densityfile,"5 30.\n");
//fprintf(densityfile,"6 \n");
}
}
if(order=="ZXY") {

for(int k=0;k<xsplit;k++)
{
  fprintf(densityfile,"#%d \n",k);
  fprintf(velocityfile,"#%d \n",k);
for(int l=0;l<ysplit;l++)
{
for(int m=0;m<zsplit;m++)
{
    fprintf(densityfile,"%d ",density[m+l*zsplit+k*ysplit*zsplit]);
    fprintf(velocityfile,"%f ",velocity[m+l*zsplit+k*ysplit*zsplit]/(density[m+l*zsplit+k*ysplit*zsplit]+0.0001));
}
fprintf(densityfile,"\n");
fprintf(velocityfile,"\n");
}
//fprintf(densityfile,"5\n");
//fprintf(densityfile,"1 10.\n");
//fprintf(densityfile,"2 15.\n");
//fprintf(densityfile,"3 20.\n");
//fprintf(densityfile,"4 25.\n");
//fprintf(densityfile,"5 30.\n");
//fprintf(densityfile,"6 \n");
}
}
if(order=="YZX") {

for(int m=0;m<zsplit;m++)
{
  fprintf(densityfile,"#%d \n",m);
  fprintf(velocityfile,"#%d \n",m);
for(int k=0;k<xsplit;k++)
{
for(int l=0;l<ysplit;l++)
{
    fprintf(densityfile,"%d ",density[m+l*zsplit+k*ysplit*zsplit]);
    fprintf(velocityfile,"%f ",velocity[m+l*zsplit+k*ysplit*zsplit]/(density[m+l*zsplit+k*ysplit*zsplit]+0.0001));
}
fprintf(densityfile,"\n");
fprintf(velocityfile,"\n");
}
//fprintf(densityfile,"5\n");
//fprintf(densityfile,"1 10.\n");
//fprintf(densityfile,"2 15.\n");
//fprintf(densityfile,"3 20.\n");
//fprintf(densityfile,"4 25.\n");
//fprintf(densityfile,"5 30.\n");
//fprintf(densityfile,"6 \n");
}  
}
fclose(densityfile);
fclose(velocityfile);



}
