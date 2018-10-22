//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################



#include "cudatools.cuh"

void compute_Zdist(void)
{
cout << type <<" is chosen"<<endl;
cout<<"Zdist: " << nElt1 <<" "<< Elt1 << " atoms are there in each structure"<<endl;

int Dirn,split;
float min,max;

if(set_max_z == "yes")
{
  minz = set_minz;
  maxz = set_maxz;
  Dirn = 2;
  min=minz;
  max=maxz;
  split=zsplit;
  
}

else if(set_max_y == "yes")
{
  minz = set_miny;
  maxz = set_maxy;
  Dirn=1;
  min=miny;
  max=maxy;
  split=ysplit;
}
else if(set_max_x == "yes")
{
  minz = set_minx;
  maxz = set_maxx;
  Dirn=0;
  min=minx;
  max=maxx;
  split=xsplit;
}



printf("%f %f\n",min,max);

    zrange=max-min;

    zint = zrange/split;

    ztick=(float*) malloc (sizeof(float)*split*2);
    
    for(int i=0;i<split;i++)
    {
      ztick[i*2]=min+(i*zint);
      ztick[i*2+1]=min+((i+1)*zint);
    }
    

FILE *xyztick = fopen("xyztick.dat","wt");
    for(int i=0;i<split;i++)
    {
      fprintf(xyztick,"%f %f \n",ztick[i*2],ztick[i*2+1]);
    }
fclose(xyztick);

    density=(int *) malloc (sizeof(int)*split);

        for(int k=0;k<split;k++)
        {
          density[k]=0;
        }

///*


    printf("Start of cuda calculation\n");
    cudaMalloc((void **)&dev_A,sizeof(float)*(nstruct/skip)*nElt1*3);
    cudaMalloc((void **)&dev_density,sizeof(int)*split);
    cudaMalloc((void **)&dev_ztick, sizeof(float)*split*2);
    cudaMemcpy(dev_A,A,sizeof(float)*(nstruct/skip)*nElt1*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_density,density,sizeof(int)*split,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_ztick,ztick,sizeof(float)*split*2,cudaMemcpyHostToDevice);
    
    cout << "Memory for storing structure data: "<<((nstruct/skip)*nElt1*3)*4/float(1000000000)<< " Gbs" << endl;
    cout << "Memory for local storage: "<<(split+(split*2))*4/float(1000000000)<< " Gbs" << endl;   
    
    dim3 dimBlock(32,32,1);
    //dim3 dimBlock(1,1,1);
    dim3 dimGrid(((nstruct/skip)+dimBlock.x-1)/dimBlock.x,(nElt1+dimBlock.y-1)/dimBlock.y,1);
    //dim3 dimGrid(nstruct,nElt1,1);
    Xdist_calc<<<dimGrid,dimBlock>>>(dev_A,dev_density,(nstruct/skip),nElt1,split,dev_ztick,Dirn);
    cudaMemcpy(density,dev_density,sizeof(int)*split,cudaMemcpyDeviceToHost);
    
    
    xvec=0;yvec=0;zvec=0;
    for(int i=0; i<(nstruct/skip);i++)
    {
        xvec+=lattice[i*6];
        yvec+=lattice[1+i*6];
        zvec+=lattice[2+i*6];
    }
    
    xvec /= (nstruct/skip);
    yvec /= (nstruct/skip);
    zvec /= (nstruct/skip);
    
    
    
    
    
    
    
    
    float cellvolume;
    
    if(cell_type=="orthorhombic")
    {
      if(set_max_z == "yes") cellvolume = xvec*yvec*zint;
      if(set_max_x == "yes") cellvolume = yvec*zvec*zint;
      if(set_max_y == "yes") cellvolume = xvec*zvec*zint;
    }
    else if(cell_type=="monoclinic")
    {
    if(set_max_z == "yes") cellvolume= xvec*xvec*zint*0.866025403; //*sin(60) for rhombohedral 
    if(set_max_x == "yes") cellvolume= yvec*zvec*zint*0.866025403;
    if(set_max_y == "yes") cellvolume= xvec*zvec*zint*0.866025403;
    }
    
    
FILE *Zdistfile=fopen("zdist.data","wt");
FILE *Zdistdenfile=fopen("zdistden.data","wt");
int zstart = 0;
int zend = split;
//int startfromzero=0;  not used because relative distances get messed up
for(int m=zstart;m<zend;m++)
{
  //if(startfromzero == 0 && float(density[m])/(xvec*yvec*zint*(nstruct/skip)) == 0) {}
  //else {fprintf(Zdistfile,"%f %f \n ",zint*startfromzero, float(density[m])/(xvec*yvec*zint*(nstruct/skip)));startfromzero++;}
  fprintf(Zdistfile,"%f %f \n ",ztick[m*2+1], float(density[m])/(cellvolume*(nstruct/skip)));
  fprintf(Zdistdenfile,"%f %f \n ",ztick[m*2+1], float(density[m]));
}
fclose(Zdistfile);
fclose(Zdistdenfile);
cudaFree(dev_density);
cudaFree(dev_A);
cudaFree(dev_ztick);

}
