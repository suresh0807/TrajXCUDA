//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"




void compute_diffusion_lifetime()
{
  
    select_atoms(msd_for);
  
/*
  cout << "metal "<< metal_species<<" counts "<<count_metal<<endl;
  
  int *exch;
  int *dev_exch;
  exch= (int *) malloc (sizeof(int)*nElt1*(nstruct/skip));
  Aintnum = (int *) malloc (sizeof(int)*(nstruct/skip));
  Abulknum = (int *) malloc (sizeof(int)*(nstruct/skip));
// Initialize distance matrix and histogram matrix
  for(int l=0; l<nstruct/skip; l++)
    {
      for(int i=0; i< nElt1; i++)
        {
          exch[i+l*nElt1]=0;
        }     
        Abulknum[l]=0;
        Aintnum[l]=0;
      }

//Allocate memory in GPU device

    cudaMalloc((void **)&dev_B,sizeof(float)*count_metal*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_A,sizeof(float)*nElt1*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_exch,sizeof(int)*nElt1*(nstruct/skip));
    cudaMalloc((void **)&dev_lattice,sizeof(float)*6*(nstruct/skip));

//Copy data from host to device

    cudaMemcpy(dev_B,METAL,sizeof(float)*count_metal*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_A,A,sizeof(float)*nElt1*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    

//dim3 struct to define elements of the execution configuration


    dim3 dimBlocke(10,10,10);
    dim3 dimGride((nElt1+dimBlocke.x-1)/dimBlocke.x,(count_metal+dimBlocke.y-1)/dimBlocke.y,((nstruct/skip)+dimBlocke.z-1)/dimBlocke.z);

//Cuda kernal execution for distance matrix with CUDA timing API commands

    //For surface metal atoms
    cudaMemcpy(dev_exch,exch,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice);
    cudaMemcpy(dev_lattice,lattice,sizeof(int)*6*(nstruct/skip),cudaMemcpyHostToDevice);
    if (cell_type == "orthorhombic") 
    {
    covermat<<<dimGride,dimBlocke>>>(dev_A,dev_B,dev_exch,bondist,nElt1,count_metal,(nstruct/skip),dev_lattice);
    }
    else if (cell_type == "monoclinic")
    {
     covermatmono<<<dimGride,dimBlocke>>>(dev_A,dev_B,dev_exch,bondist,nElt1,count_metal,(nstruct/skip),dev_lattice);
      //covermat<<<dimGride,dimBlocke>>>(dev_A,dev_B,dev_exch,bondist,nElt1,count_metal,(nstruct/skip),dev_lattice);
    }
    cudaMemcpy(exch,dev_exch,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyDeviceToHost);
    //For substrate atoms
     cudaFree(dev_A);
    cudaFree(dev_B);
    cudaFree(dev_lattice);
    //cudaFree(dev_exch);
    for(int l=0; l<nstruct/skip; l++)
    {
      for(int i=0; i< nElt1; i++)
        {
          if (exch[i+l*nElt1]>0) {exch[i+l*nElt1]=1;}
        }     
      }
    cudaMemcpy(dev_exch,exch,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice);
    ofstream exchplot;
    exchplot.open("exchange.data");

 for(int l=0; l<(nstruct/skip); l++)
    {
      exchplot <<l<<" ";
      for(int j=0; j<nElt1; j++)
      {
      
          exchplot << exch[j+l*nElt1]<<" ";
      }
      exchplot <<endl;
    }
    exchplot.close();    

 //ofstream intnum;
 //intnum.open("intnum.data");
 
      for(int l=0; l<nstruct/skip; l++)
      {
      for(int i=0; i< nElt1; i++)
        {
          if (exch[i+l*nElt1]>0) {exch[i+l*nElt1]=1; Aintnum[l]++;}
          else if (exch[i+l*nElt1]==0) Abulknum[l]++;
        }
      //for(int i=0; i< nElt2; i++)
      //  {
      //    if (exch2[i+l*nElt2]>0) {exch2[i+l*nElt2]=1; Bintnum[l]++;}
      //    else if (exch2[i+l*nElt2]==0) Bbulknum[l]++;
      //  }
      //intnum << l<<" "<<Aintnum[l]<<endl;
      }
  //intnum.close();
  */
///////////////////////unwrap after exch matrix - sensible///////////////
//##########################################################################  

if( unwrap == "yes")
{

cout << "Unwrapping the trajectory begins"<<endl;
  float lr[3];

for(int i =0; i<(nstruct/skip); i++)
{
 for(int j =0 ; j<nElt1; j++)
 {
  for(int k =0; k<3;k++)
  {
   Aint[j*3+i*nElt1*3+k]=0.0;  
  }
  if(i==0)
  {
   for(int k =0; k<3;k++)
   {
    Aint[j*3+i*nElt1*3+k]=A[j*3+i*nElt1*3+k];
   }
  }
 }
}

if(cell_type == "orthorhombic")
 {
  for(int i=1; i<(nstruct/skip); i++)
  {//cout<<lnElt1<<endl;
  //cout<<"hi"<<endl;
      for(int j=0;j<nElt1;j++)
      {//cout<<Elt1<<" ";
          for(int k=0 ;k<3;k++)
          {
          lr[k] = A[j*3+i*nElt1*3+k] - A[j*3+(i-1)*nElt1*3+k];
          if(abs(lr[k]) > lattice[k+i*6]/2.0)
          {
               //cout<<"I am working"<<endl;
              if(lr[k] > 0)
              {
              lr[k] = abs(lr[k]) - lattice[k+i*6];
              Aint[j*3+i*nElt1*3+k] = Aint[j*3+(i-1)*nElt1*3+k] - abs(lr[k]); 
              }
              else
              {
              lr[k] = abs(lr[k]) - lattice[k+i*6];
              Aint[j*3+i*nElt1*3+k] = Aint[j*3+(i-1)*nElt1*3+k] + abs(lr[k]); 
              }     
          }
           else
           {
              Aint[j*3+i*nElt1*3+k] = Aint[j*3+(i-1)*nElt1*3+k] + lr[k]; 
           } //cout<<Aint[j*3+i*lnElt1*3+k]<<" ";
          }//cout<<endl;
      }
  }
 }
 else if(cell_type == "monoclinic")
 {
  for(int i=1; i<(nstruct/skip); i++)
  {
      for(int j=0;j<nElt1;j++)
      {        
        lr[0]=A[j*3+i*nElt1*3]-A[j*3+(i-1)*nElt1*3];
        lr[1]=A[j*3+i*nElt1*3+1]-A[j*3+(i-1)*nElt1*3+1];
        lr[2]=A[j*3+i*nElt1*3+2]-A[j*3+(i-1)*nElt1*3+2];
	
	if(abs(lr[1]) > lattice[1+i*6]/2.0)
        {
        if(lr[1] > 0)
        {
            lr[1] = abs(lr[1]) - lattice[1+i*6];
            Aint[j*3+i*nElt1*3+1] = Aint[j*3+(i-1)*nElt1*3+1] - abs(lr[1]);
            lr[0] = abs(lr[0]) - lattice[3+i*6];
            Aint[j*3+i*nElt1*3+0] = Aint[j*3+(i-1)*nElt1*3+0] - abs(lr[0]);
	   // Aint[j*3+i*lnElt1*3+0] = Aint[j*3+i*lnElt1*3+0] - lattice[3+i*6];
        }
        else
        {
            lr[1] = abs(lr[1]) - lattice[1+i*6];
            Aint[j*3+i*nElt1*3+1] = Aint[j*3+(i-1)*nElt1*3+1] + abs(lr[1]);
            lr[0] = abs(lr[0]) - lattice[3+i*6];
            Aint[j*3+i*nElt1*3+0] = Aint[j*3+(i-1)*nElt1*3+0] + abs(lr[0]);
	    //Aint[j*3+i*lnElt1*3+0] = Aint[j*3+i*lnElt1*3+0] + lattice[3+i*6];
        }
        }
        else
        {
        Aint[j*3+i*nElt1*3+1] = Aint[j*3+(i-1)*nElt1*3+1] + lr[1];
        }
	
	
        for(int k=0;k<3;k=k+2)
        {
        if(abs(lr[k]) > lattice[k+i*6]/2.0)
        {
            if(lr[k] > 0)
            {
                lr[k]= abs(lr[k]) - lattice[k+i*6];
                Aint[j*3+i*nElt1*3+k] = Aint[j*3+(i-1)*nElt1*3+k] - abs(lr[k]);
            }
            else
            {
                lr[k]= abs(lr[k]) - lattice[k+i*6];
                Aint[j*3+i*nElt1*3+k] = Aint[j*3+(i-1)*nElt1*3+k] + abs(lr[k]);
            }
        }
        else
        {
            Aint[j*3+i*nElt1*3+k] = Aint[j*3+(i-1)*nElt1*3+k] + lr[k];
        }
        }
      }
  }
 }


for(int i =0; i<(nstruct/skip); i++)
{
for(int j =0 ; j<nElt1; j++)
{
for(int k =0; k<3;k++)
{
A[j*3+i*nElt1*3+k]=Aint[j*3+i*nElt1*3+k];
}
}
}

if(unwrapout=="yes")
{

ofstream unwrapped;
unwrapped.open("input_unwrapped.xyz");

for(int i =0; i<(nstruct/skip); i++)
{
unwrapped<<nElt1<<endl;
unwrapped<<endl;
for(int j =0 ; j<nElt1; j++)
{unwrapped<<Elt1 <<" ";
for(int k =0; k<3;k++)
{
unwrapped<<A[j*3+i*nElt1*3+k]<<" ";
}
unwrapped<<endl;
}
}
unwrapped.close();
}

free(Aint);
}
//#########################################################################
/////////////////////////////unwrapper finished////////////////////////////
  
  
  
    float intnum=0.0,bulknum=0.0;
    for(int l=0; l<nstruct/skip; l++)
    {
      intnum+=Aintnum[l];
      bulknum+=Abulknum[l];
    }
    intnum/=(nstruct/skip);
    bulknum/=(nstruct/skip);
  
    cout << intnum << " "<< bulknum<<endl;
    int DOF;
    if(diffuse_direction == "xyz") {xsrt=0;xend=3;xski=1;DOF=6;}
    else if(diffuse_direction == "xy") {xsrt=0;xend=2;xski=1;DOF=4;}
    else if(diffuse_direction == "x") {xsrt=0;xend=1;xski=1;DOF=2;}
    else if(diffuse_direction == "y") {xsrt=1;xend=2;xski=1;DOF=2;}
    else if(diffuse_direction == "z") {xsrt=2;xend=3;xski=1;DOF=2;}
    else if(diffuse_direction == "xz") {xsrt=0;xend=3;xski=2;DOF=4;}
    else if(diffuse_direction == "yz") {xsrt=1;xend=3;xski=1;DOF=4;}
 

    total_time = (nstruct)*timestep;
    cout << "Each frame resolves to "<<timestep<< " fs evolution"<<endl;
    cout << "Total simulation time : " << total_time/1000000 <<" ns"<< endl;
    num_bins =  int ((nstruct )/((diffuse_time*1000)/timestep));
    num_bins -= 1; 
    cout << "The analysis will be done for every "<< diffuse_time <<" ps"<<endl;
    cout <<num_bins <<" number of bins are chosen for the MSD analysis"<<endl;
    SD_store = int((diffuse_time*1000)/timestep);
    origins = SD_store;
    cout <<"This will correlate for "<<SD_store<<" frames in the input trajectory "<<endl;
    cout <<origins/skips<<" restarts from each bin will be taken"<<endl;


    SD=(float *) malloc (sizeof(float)*nElt1*SD_store*origins/skips);
    for(int i=0;i<origins/skips;i++)
    {
      for(int j=0;j<SD_store;j++)
      {
	for(int k=0; k<nElt1; k++)
	{
          SD[k+j*nElt1+i*nElt1*SD_store]=0.0;
	}
      }
    }
        SDsum=(float *) malloc (sizeof(float)*nElt1*SD_store*origins/skips);
    for(int i=0;i<origins/skips;i++)
    {
      for(int j=0;j<SD_store;j++)
      {
	for(int k=0; k<nElt1; k++)
	{
          SDsum[k+j*nElt1+i*nElt1*SD_store]=0.0;
	}
      }
    }
///*

    
    
    printf("Start of cuda calculation\n");
    cudaMalloc((void **)&dev_A,sizeof(float)*nstruct*nElt1*3);
    cudaMalloc((void **)&dev_B,sizeof(float)*nstruct*nElt2*3);
    cudaMalloc((void **)&dev_SD,sizeof(float)*nElt1*SD_store*origins/skips);
    cudaMemcpy(dev_A,A,sizeof(float)*(nstruct/skip)*nElt1*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_B,B,sizeof(float)*(nstruct/skip)*nElt2*3,cudaMemcpyHostToDevice);
    cout << "Memory required to store coordinate information: "<<nstruct*nElt1*3*2*4/float(1000000000)<< " Gbs " << endl;
    cout << "Memory for local MSD storage: "<<((nElt1*nstruct/skips)+(nElt1*SD_store))*4/float(1000000000)<< " Gbs" << endl;


    
   // cout << "Memory required to store coordinate information: "<<nstruct*nElt1*3*2*4/float(1000000000)<< " Gbs " << endl;
   // cout << "Memory for local MSD storage: "<<nElt1*SD_store*origins/skips*4/float(1000000000)<< " Gbs" << endl;


    dim3 dimBlock(32,1,32);
    dim3 dimGrid((nElt1+dimBlock.x-1)/dimBlock.x,1,((origins/skips)+dimBlock.z-1)/dimBlock.z);

    //dim3 dimBlock(10,10,10);
    //dim3 dimGrid((nElt1+dimBlock.x-1)/dimBlock.x,(SD_store+dimBlock.y-1)/dimBlock.y,((origins/skips)+dimBlock.z-1)/dimBlock.z);    
    
  //    cudaMemcpy(dev_SD,SD,sizeof(float)*nElt1*SD_store*origins,cudaMemcpyHostToDevice);    

int *numc, *dev_numc;
numc = (int *) malloc (sizeof(int)*SD_store*origins/skips);
cudaMalloc((void **)&dev_numc,sizeof(int)*SD_store*origins/skips);
    
int *inout, *dev_inout;
inout = (int *) malloc (sizeof(int)*nElt1*origins/skips);
cudaMalloc((void **)&dev_inout,sizeof(int)*nElt1*origins/skips);

    for(int i=0; i<num_bins;i++)
    //for(int i=0; i<num_bins;i++)
    {
             for(int i1=0;i1<origins/skips;i1++)
              {
               for(int j1=0;j1<SD_store;j1++)
                {
		  numc[j1+i1*SD_store]=0;
	         for(int k1=0; k1<nElt1; k1++)
	          {
                    SD[k1+j1*nElt1+i1*nElt1*SD_store]=0.0;
		    inout[k1+i1*nElt1]=1;
	          }
                }
              }
   
      	      int whichwater;
      if(msd_for=="int"||msd_for=="box"||msd_for=="bulk"){whichwater=1;}
      else if(msd_for=="all"){whichwater=2;}
cudaMemcpy(dev_SD,SD,sizeof(float)*nElt1*SD_store*origins/skips,cudaMemcpyHostToDevice);
cudaMemcpy(dev_numc,numc,sizeof(int)*SD_store*origins/skips,cudaMemcpyHostToDevice);
cudaMemcpy(dev_inout,inout, sizeof(int)*nElt1*(origins/skips),cudaMemcpyHostToDevice);
int chase=1;
      for(int j=1; j<SD_store; j++)
      {
      //cudaMemcpy(dev_inout,inout, sizeof(int)*nElt1*(origins/skips),cudaMemcpyHostToDevice);
      SD_calc<<<dimGrid,dimBlock>>>(dev_A,dev_B,dev_SD,dev_exch,dev_numc,dev_inout,SD_store,nElt1,i,j,origins,skips,xsrt,xend,xski,whichwater,chase);
      //SD_calc<<<dimGrid,dimBlock>>>(dev_A,dev_B,dev_SD,dev_exch,dev_numc,dev_inout,SD_store,nElt1,i,j,origins,skips,xsrt,xend,xski,whichwater);
      //cudaMemcpy(inout,dev_inout, sizeof(int)*nElt1*(origins/skips),cudaMemcpyHostToDevice);
      }
cudaMemcpy(SD,dev_SD,sizeof(float)*nElt1*SD_store*origins/skips,cudaMemcpyDeviceToHost);
cudaMemcpy(numc,dev_numc,sizeof(int)*SD_store*origins/skips,cudaMemcpyDeviceToHost);

    for(int i1=0;i1<origins/skips;i1++)
    {
     for(int j1=1;j1<SD_store;j1++)
     {
      for(int k1=0; k1<nElt1; k1++)
      {
       SD[k1+j1*nElt1+i1*nElt1*SD_store]/=float(numc[j1+i1*SD_store]);
      }
     }
    }
    
    for(int i1=0;i1<origins/skips;i1++)
    {
      for(int j1=0;j1<SD_store;j1++)
      {
	for(int k1=0; k1<nElt1; k1++)
	{
          SDsum[k1+j1*nElt1+i1*nElt1*SD_store]+=SD[k1+j1*nElt1+i1*nElt1*SD_store];
	}
      }
    }
     cout << (i+1)*diffuse_time << " pico seconds done"<<endl;  
 }
    
//      cudaMemcpy(SD,dev_SD,sizeof(float)*nElt1*SD_store*origins,cudaMemcpyDeviceToHost);



    cudaFree(dev_A);
    cudaFree(dev_exch);
    cudaFree(dev_B);
    cudaFree(dev_SD);
    cudaFree(dev_numc);
    cudaFree(dev_inout);
     
    SDavgf=(float *) malloc (sizeof(float)*SD_store*nElt1);
    for(int i=0; i<SD_store; i++)
    {
      for(int j =0; j<nElt1;j++)
      {
        SDavgf[j+i*nElt1]=0.0;
      }
    }
    cudaMalloc((void **)&dev_SDsum,sizeof(float)*nElt1*SD_store*origins/skips);
    cudaMalloc((void **)&dev_SDavgf,sizeof(float)*SD_store*nElt1);
    cudaMemcpy(dev_SDavgf,SDavgf,sizeof(float)*SD_store*nElt1,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_SDsum,SDsum,sizeof(float)*nElt1*SD_store*origins/skips,cudaMemcpyHostToDevice);
    
    dim3 dimBlockf(32,32,1);
    dim3 dimGridf((nElt1+dimBlockf.x-1)/dimBlockf.x,(SD_store+dimBlockf.y-1)/dimBlockf.y,1);
    cout<<"start reduce "<< origins/skips <<" "<<nElt1<< " "<<SD_store<<endl;
    
    SDreducef<<<dimGridf,dimBlockf>>>(dev_SDsum,dev_SDavgf,nElt1,SD_store,origins,skips);
    
    cudaMemcpy(SDavgf,dev_SDavgf,sizeof(float)*SD_store*nElt1,cudaMemcpyDeviceToHost);
    cudaFree(dev_SDsum);

    int fairy=0;
    SDavg=(float *) malloc (sizeof(float)*SD_store);
    for(int j =0; j<SD_store;j++)
    {
      SDavg[j]=0.0;
    }
    cudaMalloc((void **)&dev_SDavg,sizeof(float)*SD_store);
    cudaMemcpy(dev_SDavg,SDavg,sizeof(float)*SD_store,cudaMemcpyHostToDevice);

    dim3 dimBlocks(1024,1,1);
    dim3 dimGrids((SD_store+dimBlocks.x-1)/dimBlocks.x,1,1);

    SDreduce<<<dimGrids,dimBlocks>>>(dev_SDavgf,dev_SDavg,SD_store,nElt1,fairy);
    cudaMemcpy(SDavg,dev_SDavg,sizeof(float)*SD_store,cudaMemcpyDeviceToHost);
    
    cudaFree(dev_SDavgf);
    cudaFree(dev_SDavg);
    
 /*
    if(msd_for == "int")
    {
     for(int i =1; i<SD_store;i++)
    {
      SDavg[i]=(SDavg[i]*nElt1)/ intnum;
    }
    }
    else if(msd_for == "bulk")
    {
     for(int i =1; i<SD_store;i++)
    {
      SDavg[i]=(SDavg[i]*nElt1)/ bulknum;
    }
    }      
   */ 
// printing the msd data to be visualized

    FILE *MSDplot=fopen("msd.data","wt");
      fprintf(MSDplot,"# Time (ps) MSD (A^(2)) \n");
    for(int i =1; i<SD_store;i++)
    {
      fprintf(MSDplot,"%f %f\n",(float(i)*timestep)/1000.0,SDavg[i]/num_bins);
    }
    fclose(MSDplot);
   
    FILE *DIFFplot=fopen("diffco.data","wt");
      fprintf(DIFFplot,"# Time (ps) D (10^-5 cm^2/s) \n");
    for(int i =1; i<SD_store;i++)
    {
      fprintf(DIFFplot,"%f %f\n",(float(i)*timestep)/1000.0,(SDavg[i]/(num_bins*DOF*i*timestep))*10000);
    }
    fclose(DIFFplot);


 
}
