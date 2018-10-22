//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################


#include "cudatools.cuh"

void compute_tetrawater(void)
{


  float *Onum; //tetrahedral parameter for each O atom
  int *Hbondneigh; // Ids of the O neighbours 
  float *Hbondneighdist; // OO distance for all neighbours.
  int *Oneigh;// Ids of 4 O neighbours
  float *Hdistmat, *dev_Hdistmat;// Distance matrix to capture the O within x \AA of O atom.
  
  int OHsamples=20; /// maximum number of OH bonds considered within bonddist sphere
  
    Onum = (float *) malloc (sizeof(float)*nElt1*(nstruct/skip));
    Hbondneigh = (int *) malloc (sizeof(int)*OHsamples*nElt1*(nstruct/skip));
    Oneigh = (int *) malloc (sizeof(int)*4*nElt1*(nstruct/skip));
    Hbondneighdist = (float *) malloc (sizeof(float)*OHsamples*nElt1*(nstruct/skip));
    Hdistmat = (float *) malloc (sizeof(float)*nElt1*nElt1);
    
  for(int i=0;i<nstruct/skip;i++)
  {
    for(int j=0; j<nElt1; j++)
    {
      Onum[j+i*nElt1] =0;
      for(int k =0 ;k<OHsamples;k++)
      {
      Hbondneigh[j*OHsamples+i*nElt1*OHsamples+k]=0;
      Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k]=0.0;
      }
    }
  }   
      
 for(int i=0;i<nElt1;i++)
  {
    for(int j=0; j<nElt1; j++)
    {
      Hdistmat[j+i*nElt1]=0.0;
    }
  }
      
cout <<"gpu begins distance computation"<<endl;
//Allocate memory in GPU device

    cudaMalloc((void **)&dev_A,sizeof(float)*nElt1*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_B,sizeof(float)*nElt1*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_Hdistmat,sizeof(float)*nElt1*nElt1);
    cudaMalloc((void **)&dev_lattice,sizeof(float)*6*(nstruct/skip));
//Copy data from host to device

    cudaMemcpy(dev_A,A,sizeof(float)*nElt1*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_B,A,sizeof(float)*nElt1*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_Hdistmat,Hdistmat,sizeof(float)*nElt1*nElt1,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_lattice,lattice,sizeof(int)*6*(nstruct/skip),cudaMemcpyHostToDevice);

   //dim3 struct to define elements of the execution configuration


    dim3 dimBlock(32,32,1);
    dim3 dimGrid((nElt1+dimBlock.x-1)/dimBlock.x,(nElt1+dimBlock.y-1)/dimBlock.y,1);
for(int i = 0; i<(nstruct/skip); i++) /// go into each frame
{
  
   for(int i1=0;i1<nElt1;i1++)///initialize Hdistmat for 1 frame
   {
    for(int j1=0; j1<nElt1; j1++)
    {
      Hdistmat[j1+i1*nElt1]=0.0;
    }
   }
  
    cudaMemcpy(dev_Hdistmat,Hdistmat,sizeof(float)*nElt1*nElt1,cudaMemcpyHostToDevice);
//Cuda kernal execution for distance matrix with CUDA timing API commands
    if(cell_type == "orthorhombic")
    {
          Hbondmat<<<dimGrid,dimBlock>>>(dev_A,dev_B,dev_Hdistmat,bondist,nElt1,nElt1,dev_lattice,i);
    }
    else if(cell_type == "monoclinic")
    {
          Hbondmatmono<<<dimGrid,dimBlock>>>(dev_A,dev_B,dev_Hdistmat,bondist,nElt1,nElt1,dev_lattice,i);      
    }

cudaMemcpy(Hdistmat,dev_Hdistmat,sizeof(float)*nElt1*nElt1,cudaMemcpyDeviceToHost);
   
  for(int i1=0;i1<nElt1;i1++) //oxygen1
  {
    int k=0;
    for(int j1=0; j1<nElt1; j1++) //oxygen2 
    { 
      if(Hdistmat[j1+i1*nElt1] !=0.0) // check if an oxygen atom is within the 'bonddist' distance from another O atom
      {
	Onum[i1+i*nElt1]++;
	Hbondneigh[i1*OHsamples+i*nElt1*OHsamples+k] = j1;
	Hbondneighdist[i1*OHsamples+i*nElt1*OHsamples+k] = Hdistmat[j1+i1*nElt1];
	k++;
      }
    }
  }

}   
    
  
    cudaFree(dev_A);
    cudaFree(dev_B);
    cudaFree(dev_Hdistmat);
    cudaFree(dev_lattice);


//Hbondneigh contain the IDs of 20 possible oxygen neighbour atoms.
//Onum has the number of oxygen atoms within a specified distance (see input file).

cout <<"gpu ends"<<endl; 

      
//swapping the indices of nearest oxygens to get the two shortest OH bonds up front, 
//since it usually describes the H bond

cout<<"Swapping of OO bonds based on distance "<<endl;

float swapa,swapb;
    for(int i=0;i<(nstruct/skip);i++)
     {
      for(int j=0; j<nElt1; j++)
       {
	 //if(i==0){cout<<j<<" ";}
	for(int swf=0; swf<Onum[j+i*nElt1]-1;swf++)
	{
	for(int sws=swf+1; sws<Onum[j+i*nElt1];sws++)
	{
	if(Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+swf] > Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+sws] && sws !=swf && Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+sws] !=0.0 && Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+swf] !=0.0) 
	{
	   swapa = Hbondneigh[j*OHsamples+i*nElt1*OHsamples+swf]; 
	    Hbondneigh[j*OHsamples+i*nElt1*OHsamples+swf]=Hbondneigh[j*OHsamples+i*nElt1*OHsamples+sws];
	    Hbondneigh[j*OHsamples+i*nElt1*OHsamples+sws]=swapa;
	    swapb = Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+swf]; 
	    Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+swf]=Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+sws];
	    Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+sws]=swapb;
	}
	}
	}
	for(int k =0 ;k<20;k++)
        {
	//if(i==0){cout<<Hbondneigh[j*OHsamples+i*nElt1*OHsamples+k]<<" "<<Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k]<<"  ";}
	}
	//if(i==0){cout<<endl;}
       }
       
     }	
   
cout<<" done !!!"<<endl;
   
   
  free(Hbondneighdist);
  free(Hdistmat);
  free(Onum);
  
  for(int i=0;i<nstruct/skip;i++)
  {
    for(int j=0; j<nElt1; j++)
    {
      for(int k =0 ;k<4;k++)
      {
      Oneigh[j*4+i*nElt1*4+k]=Hbondneigh[j*OHsamples+i*nElt1*OHsamples+k];
      }
    }
  }
  
  
  
  free(Hbondneigh);
  
  float angle_final;
  float *OneighXYZ, *angXYZ, *Otetra;
  OneighXYZ = (float *) malloc (sizeof(float)*5*3);
  Otetra = (float *) malloc (sizeof(float)*(nstruct/skip)*nElt1);
  angXYZ = (float *) malloc (sizeof(float)*6);
  for (int i=0;i<15;i++)
  {
    OneighXYZ[i]=0;
  }
  for (int i=0;i<6;i++)
  {
    angXYZ[i]=0;
  }
  
  for(int i=0;i<nstruct/skip;i++)
  {
    for(int j=0; j<nElt1; j++)
    {
      angle_final=0.0;
      Otetra[j+i*nElt1]=0;
      for(int ss=0; ss<3;ss++)
      {
	OneighXYZ[ss]=A[j*3+i*nElt1*3+ss];
      }
      for(int k =0 ;k<4;k++)
      {
	for(int ss=0; ss<3;ss++)
        {
	  OneighXYZ[3+k*3+ss]=A[Oneigh[j*4+i*nElt1*4+k]*3+i*nElt1*3+ss];
	}
      }
      angXYZ[0]=angle(OneighXYZ[0],OneighXYZ[1],OneighXYZ[2],OneighXYZ[3],OneighXYZ[4],OneighXYZ[5],\
	OneighXYZ[6],OneighXYZ[7],OneighXYZ[8]);
      angXYZ[1]=angle(OneighXYZ[0],OneighXYZ[1],OneighXYZ[2],OneighXYZ[3],OneighXYZ[4],OneighXYZ[5],\
	OneighXYZ[9],OneighXYZ[10],OneighXYZ[11]);
      angXYZ[2]=angle(OneighXYZ[0],OneighXYZ[1],OneighXYZ[2],OneighXYZ[3],OneighXYZ[4],OneighXYZ[5],\
	OneighXYZ[12],OneighXYZ[13],OneighXYZ[14]);
      angXYZ[3]=angle(OneighXYZ[0],OneighXYZ[1],OneighXYZ[2],OneighXYZ[6],OneighXYZ[7],OneighXYZ[8],\
	OneighXYZ[9],OneighXYZ[10],OneighXYZ[11]);
      angXYZ[4]=angle(OneighXYZ[0],OneighXYZ[1],OneighXYZ[2],OneighXYZ[6],OneighXYZ[7],OneighXYZ[8],\
	OneighXYZ[12],OneighXYZ[13],OneighXYZ[14]);
      angXYZ[5]=angle(OneighXYZ[0],OneighXYZ[1],OneighXYZ[2],OneighXYZ[9],OneighXYZ[10],OneighXYZ[11],\
	OneighXYZ[12],OneighXYZ[13],OneighXYZ[14]);
      for(int ss=0;ss<6;ss++)
      {
      angle_final= angle_final + square(cos(angXYZ[ss])+(1.0/3.0));
      }
      Otetra[j+i*nElt1]= 1.0- ( (3.0/8.0) * angle_final );
      //if(i==0 && j==0) {cout<<angXYZ[0]<<" "<<cos(angXYZ[0])<<" "<<angle_final<<" "<<Otetra[j+i*nElt1]<<endl;}
    }
  }
  
 /* 
  for(int i=0;i<1;i++)
  {
    for(int j=0;j<nElt1;j++)
    {
      cout<<j<<" "<<Oneigh[j*4+i*nElt1*4+0]<<" "<<Oneigh[j*4+i*nElt1*4+1]<<" "<<\
      Oneigh[j*4+i*nElt1*4+2]<<" "<<Oneigh[j*4+i*nElt1*4+3]<<" "<<Otetra[j+i*nElt1]<<endl;
    }
  }
  */
 
 
  free(OneighXYZ);
  free(Oneigh);
  free(angXYZ);
  

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

float *tetradensity;
    density=(int *) malloc (sizeof(int)*split);
tetradensity=(float *) malloc (sizeof(float)*split);
        for(int k=0;k<split;k++)
        {
          density[k]=0;
	  tetradensity[k]=0;
        }

///*

    float *dev_tetra, *dev_tetradensity;
    
    printf("Start of cuda calculation\n");
    cudaMalloc((void **)&dev_A,sizeof(float)*(nstruct/skip)*nElt1*3);
    cudaMalloc((void **)&dev_tetra,sizeof(float)*(nstruct/skip)*nElt1);
    cudaMalloc((void **)&dev_density,sizeof(int)*split);
    cudaMalloc((void **)&dev_tetradensity,sizeof(float)*split);
    cudaMalloc((void **)&dev_ztick, sizeof(float)*split*2);
    cudaMemcpy(dev_A,A,sizeof(float)*(nstruct/skip)*nElt1*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_tetra,Otetra,sizeof(float)*(nstruct/skip)*nElt1,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_density,density,sizeof(int)*split,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_tetradensity,tetradensity,sizeof(float)*split,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_ztick,ztick,sizeof(float)*split*2,cudaMemcpyHostToDevice);
    
    cout << "Memory for storing structure data: "<<((nstruct/skip)*nElt1*3)*4/float(1000000000)<< " Gbs" << endl;
    cout << "Memory for local storage: "<<(split+(split*2))*4/float(1000000000)<< " Gbs" << endl;   
    
    dim3 dimBlockr(32,32,1);
    //dim3 dimBlock(1,1,1);
    dim3 dimGridr(((nstruct/skip)+dimBlock.x-1)/dimBlock.x,(nElt1+dimBlock.y-1)/dimBlock.y,1);
    //dim3 dimGrid(nstruct,nElt1,1);
    Xdist_calc<<<dimGridr,dimBlockr>>>(dev_A,dev_tetra,dev_tetradensity,dev_density,(nstruct/skip),nElt1,split,dev_ztick,Dirn);
    cudaMemcpy(density,dev_density,sizeof(int)*split,cudaMemcpyDeviceToHost);
    cudaMemcpy(tetradensity,dev_tetradensity,sizeof(float)*split,cudaMemcpyDeviceToHost);
    
    
    
FILE *Zdistfile=fopen("tetradist.data","wt");
int zstart = 0;
int zend = split;
//int startfromzero=0;  not used because relative distances get messed up
for(int m=zstart;m<zend;m++)
{
  //if(startfromzero == 0 && float(density[m])/(xvec*yvec*zint*(nstruct/skip)) == 0) {}
  //else {fprintf(Zdistfile,"%f %f \n ",zint*startfromzero, float(density[m])/(xvec*yvec*zint*(nstruct/skip)));startfromzero++;}
  fprintf(Zdistfile,"%f %f \n ",ztick[m*2+1], float(tetradensity[m]/density[m]));
}
fclose(Zdistfile);
cudaFree(dev_density);
cudaFree(dev_tetradensity);
cudaFree(dev_A);
cudaFree(dev_tetra);
cudaFree(dev_ztick);
  
   free(Otetra);
   free(tetradensity);
}
