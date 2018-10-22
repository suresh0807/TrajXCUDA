//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################



#include "cudatools.cuh"

void compute_Vdist(void)
{
cout << type <<" is chosen"<<endl;
cout<<"Vdist: " << nElt1 <<" "<< Elt1 << " atoms are there in each structure"<<endl;

int Dirn;
float swp;

if(set_max_z == "yes")
{
  minz = set_minz;
  maxz = set_maxz;
  Dirn = 2;
}

else if(set_max_y == "yes")
{
  minz = set_miny;
  maxz = set_maxy;
  Dirn=1;
  swp=zvec;
  zvec=yvec;
  yvec=swp;
}
else if(set_max_x == "yes")
{
  minz = set_minx;
  maxz = set_maxx;
  Dirn=0;
  swp=zvec;
  zvec=xvec;
  xvec=swp;
}
printf("%f %f\n",minz,maxz);

    zrange=maxz-minz;

    zint = zrange/zsplit;

    ztick=(float*) malloc (sizeof(float)*zsplit*2);
    
    for(int i=0;i<zsplit;i++)
    {
      ztick[i*2]=minz+(i*zint);
      ztick[i*2+1]=minz+((i+1)*zint);
    }
    
float *Vxdensity, *Vydensity, *Vzdensity;
float *dev_Vxdensity, *dev_Vydensity, *dev_Vzdensity;
    
FILE *xyztick = fopen("xyztick.dat","wt");
    for(int i=0;i<zsplit;i++)
    {
      fprintf(xyztick,"%f %f \n",ztick[i*2],ztick[i*2+1]);
    }
fclose(xyztick);

    Vdensity=(float *) malloc (sizeof(float)*zsplit);
    Vxdensity=(float *) malloc (sizeof(float)*zsplit);
    Vydensity=(float *) malloc (sizeof(float)*zsplit);
    Vzdensity=(float *) malloc (sizeof(float)*zsplit);

        for(int k=0;k<zsplit;k++)
        {
          Vdensity[k]=0.0;
	  Vxdensity[k]=0.0;
	  Vydensity[k]=0.0;
	  Vzdensity[k]=0.0;
        }
    density=(int *) malloc (sizeof(int)*zsplit);

        for(int k=0;k<zsplit;k++)
        {
          density[k]=0;
        }
///*


    printf("Start of cuda calculation\n");
    cudaMalloc((void **)&dev_A,sizeof(float)*(nstruct/skip)*nElt1*3);
    cudaMalloc((void **)&dev_VEL,sizeof(float)*(nstruct/skip)*nElt1*3);
    cudaMalloc((void **)&dev_Vdensity,sizeof(float)*zsplit);
    cudaMalloc((void **)&dev_Vxdensity,sizeof(float)*zsplit);
    cudaMalloc((void **)&dev_Vydensity,sizeof(float)*zsplit);
    cudaMalloc((void **)&dev_Vzdensity,sizeof(float)*zsplit);
    cudaMalloc((void **)&dev_ztick, sizeof(float)*zsplit*2);
    cudaMemcpy(dev_A,A,sizeof(float)*(nstruct/skip)*nElt1*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_VEL,VEL,sizeof(float)*(nstruct/skip)*nElt1*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_Vdensity,Vdensity,sizeof(float)*zsplit,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_Vxdensity,Vxdensity,sizeof(float)*zsplit,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_Vydensity,Vydensity,sizeof(float)*zsplit,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_Vzdensity,Vzdensity,sizeof(float)*zsplit,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_ztick,ztick,sizeof(float)*zsplit*2,cudaMemcpyHostToDevice);
    
    cout << "Memory for storing structure data: "<<((nstruct/skip)*nElt1*3)*4/float(1000000000)<< " Gbs" << endl;
    cout << "Memory for local storage: "<<((zsplit*4)+(zsplit*2))*4/float(1000000000)<< " Gbs" << endl;   
    
    //int dirn;
    
    //if(vel_dirn=="all"){dirn=0;}
    //else if(vel_dirn=="x"){dirn=1;}
    //else if(vel_dirn=="y"){dirn=2;}
    //else if(vel_dirn=="z"){dirn=3;}
    
    dim3 dimBlock(32,32,1);
    //dim3 dimBlock(1,1,1);
    dim3 dimGrid(((nstruct/skip)+dimBlock.x-1)/dimBlock.x,(nElt1+dimBlock.y-1)/dimBlock.y,1);
    //dim3 dimGrid(nstruct,nElt1,1);
    Vdist_calc<<<dimGrid,dimBlock>>>(dev_A,dev_VEL,dev_Vdensity,dev_Vxdensity,dev_Vydensity,dev_Vzdensity,(nstruct/skip),nElt1,zsplit,dev_ztick,Dirn);
    cudaMemcpy(Vdensity,dev_Vdensity,sizeof(float)*zsplit,cudaMemcpyDeviceToHost);
    cudaMemcpy(Vxdensity,dev_Vxdensity,sizeof(float)*zsplit,cudaMemcpyDeviceToHost);
    cudaMemcpy(Vydensity,dev_Vydensity,sizeof(float)*zsplit,cudaMemcpyDeviceToHost);
    cudaMemcpy(Vzdensity,dev_Vzdensity,sizeof(float)*zsplit,cudaMemcpyDeviceToHost);
    cudaFree(dev_VEL);
    
    printf("Start of cuda calculation\n");
    cudaMalloc((void **)&dev_density,sizeof(int)*zsplit);
    cudaMemcpy(dev_density,density,sizeof(int)*zsplit,cudaMemcpyHostToDevice);
    
    cout << "Memory for storing structure data: "<<((nstruct/skip)*nElt1*3)*4/float(1000000000)<< " Gbs" << endl;
    cout << "Memory for local storage: "<<(zsplit+(zsplit*2))*4/float(1000000000)<< " Gbs" << endl;   
    
    dim3 dimBlocks(32,32,1);
    //dim3 dimBlock(1,1,1);
    dim3 dimGrids(((nstruct/skip)+dimBlocks.x-1)/dimBlocks.x,(nElt1+dimBlocks.y-1)/dimBlocks.y,1);
    //dim3 dimGrid(nstruct,nElt1,1);
    Xdist_calc<<<dimGrids,dimBlocks>>>(dev_A,dev_density,(nstruct/skip),nElt1,zsplit,dev_ztick,Dirn);
    cudaMemcpy(density,dev_density,sizeof(int)*zsplit,cudaMemcpyDeviceToHost);
    
    
    
FILE *Vdistfile=fopen("vdist.data","wt");
//FILE *Vdistnormfile=fopen("vdistnorm.data","wt");
//FILE *Zdistfile=fopen("zdist.data","wt");
//FILE *Zdist1file=fopen("zdist1.data","wt");
int zstart = 0;
int zend = zsplit;
//int startfromzero=0;  not used because relative distances get messed up
for(int m=zstart;m<zend;m++)
{
  //if(startfromzero == 0 && float(density[m])/(xvec*yvec*zint*(nstruct/skip)) == 0) {}
  //else {fprintf(Zdistfile,"%f %f \n ",zint*startfromzero, float(density[m])/(xvec*yvec*zint*(nstruct/skip)));startfromzero++;}
  //fprintf(Zdistfile,"%f %f \n ",ztick[m*2+1], float(density[m])/(xvec*yvec*zint*(nstruct/skip)));
  //fprintf(Zdist1file,"%f %f \n ",ztick[m*2+1], float(density[m]));
  fprintf(Vdistfile,"%f %f %f %f %f \n ",ztick[m*2+1], Vdensity[m]/float(density[m]),Vxdensity[m]/float(density[m]),Vydensity[m]/float(density[m]),Vzdensity[m]/float(density[m]));
  //fprintf(Vdistnormfile,"%f %f \n ",ztick[m*2+1], float(Vdensity[m])/float(density[m]));
//  fprintf(Zdistfile,"%f %f \n ",ztick[m*2+1], float(density[m])/((nstruct/skip)));
}
//fclose(Zdistfile);
//fclose(Zdist1file);
fclose(Vdistfile);
//fclose(Vdistnormfile);

cudaFree(dev_Vdensity);
cudaFree(dev_Vxdensity);
cudaFree(dev_Vydensity);
cudaFree(dev_Vzdensity);
cudaFree(dev_density);
cudaFree(dev_A);
cudaFree(dev_ztick);

}
