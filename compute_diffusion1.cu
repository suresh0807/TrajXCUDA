//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"



void compute_diffusion1()
{
    if(plot=="int")        { nElt1 = nA_int;}
    else if(plot=="bulk")  { nElt1 = nA_bulk;}
    
    
    if(diffuse_direction == "xyz") {xsrt=0;xend=3;xski=1;}
    else if(diffuse_direction == "xy") {xsrt=0;xend=2;xski=1;}
    else if(diffuse_direction == "x") {xsrt=0;xend=1;xski=1;}
    else if(diffuse_direction == "y") {xsrt=1;xend=2;xski=1;}
    else if(diffuse_direction == "z") {xsrt=2;xend=3;xski=1;}
    else if(diffuse_direction == "xz") {xsrt=0;xend=3;xski=2;}
    else if(diffuse_direction == "yz") {xsrt=1;xend=3;xski=1;}
    
    
    total_time = (nstruct)*timestep;
    cout << "Each frame resolves to "<<timestep<< " fs evolution"<<endl;
    cout << "Total simulation time : " << total_time/1000000 <<" ns"<< endl;
    cout << "The analysis will be done for every "<< diffuse_time <<" ps"<<endl;
    SD_store = int((diffuse_time*1000)/timestep);
    cout <<"In each bin "<<SD_store<<" frames are kept"<<endl;
    cout <<nstruct/skips<<" restarts will be taken"<<endl;

    SD=(float *) malloc (sizeof(float)*nElt1*SD_store);

      for(int j=0;j<SD_store;j++)
      {
	for(int k=0; k<nElt1; k++)
	{
          SD[k+j*nElt1]=0.0;
	}
      }
      cout<<"bull"<<endl;
    SDsum=(float *) malloc (sizeof(float)*nElt1*(nstruct/skips));

      for(int j=0;j<nstruct/skips;j++)
      {
	for(int k=0; k<nElt1; k++)
	{
          SDsum[k+j*nElt1]=0.0;
	}
      }
    cout<<"bull"<<endl;
    printf("Start of cuda calculation\n");
    cudaMalloc((void **)&dev_A,sizeof(float)*(nstruct/skip)*nElt1*3);
    cudaMalloc((void **)&dev_B,sizeof(float)*(nstruct/skip)*nElt2*3);
    cudaMalloc((void **)&dev_SDsum,sizeof(float)*nElt1*nstruct/skips);
    cudaMalloc((void **)&dev_SD,sizeof(float)*nElt1*SD_store);
   
    if(plot=="int")
    {
      cudaMemcpy(dev_A,A_int,sizeof(float)*nstruct*nA_int*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_B,B_int,sizeof(float)*nstruct*nB_int*3,cudaMemcpyHostToDevice);
        cout << "Memory required to store coordinate information: "<<nstruct*nElt1*3*2*4/float(1000000000)<< " Gbs " << endl;
    cout << "Memory for local MSD storage: "<<((nElt1*nstruct/skips)+(nElt1*SD_store))*4/float(1000000000)<< " Gbs" << endl;
    }
    else if(plot=="bulk")
    {
      cudaMemcpy(dev_A,A_bulk,sizeof(float)*nstruct*nA_bulk*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_B,B_bulk,sizeof(float)*nstruct*nB_bulk*3,cudaMemcpyHostToDevice);
        cout << "Memory required to store coordinate information: "<<nstruct*nElt1*3*2*4/float(1000000000)<< " Gbs " << endl;
    cout << "Memory for local MSD storage: "<<((nElt1*nstruct/skips)+(nElt1*SD_store))*4/float(1000000000)<< " Gbs" << endl;
    }
    else if(plot=="all")
    {
cout<<"hi "<< (nstruct/skip)*nElt1*3<<endl;
      cudaMemcpy(dev_A,A,sizeof(float)*(nstruct/skip)*nElt1*3,cudaMemcpyHostToDevice);
cout<<"hi"<<endl;
    cudaMemcpy(dev_B,B,sizeof(float)*(nstruct/skip)*nElt2*3,cudaMemcpyHostToDevice);
        cout << "Memory required to store coordinate information: "<<nstruct*nElt1*3*2*4/float(1000000000)<< " Gbs " << endl;
    cout << "Memory for local MSD storage: "<<((nElt1*nstruct/skips)+(nElt1*SD_store))*4/float(1000000000)<< " Gbs" << endl;
    }
    
    


    cudaMemcpy(dev_SD,SD,sizeof(float)*nElt1*SD_store,cudaMemcpyHostToDevice);

    for(int i=1; i<SD_store;i++)
    {      
      
      for(int j1=0;j1<nstruct/skips;j1++)
      {
	for(int k1=0; k1<nElt1; k1++)
	{
          SDsum[k1+j1*nElt1]=0.0;
	}
      }
     
      cudaMemcpy(dev_SDsum,SDsum,sizeof(float)*nElt1*nstruct/skips,cudaMemcpyHostToDevice);

      //cout<<"Ha start"<<endl;

      dim3 dimBlock(32,1,32);
      dim3 dimGrid((nElt1+dimBlock.x-1)/dimBlock.x,1,(nstruct+dimBlock.z-1)/dimBlock.z);
      SD_calc<<<dimGrid,dimBlock>>>(dev_A,dev_B,dev_SDsum,SD_store,nElt1,i,(nstruct-i)/skips,skips,xsrt,xend,xski);

      //cout<<"Ha end"<<endl;
      
      dim3 dimBlockf(1024,1,1);
      dim3 dimGridf((nElt1+dimBlockf.x-1)/dimBlockf.x,1,1);      
      SDreducef<<<dimGridf,dimBlockf>>>(dev_SDsum,dev_SD,nElt1,SD_store,(nstruct-i)/skips,skips,i);
      
      
      cudaMemcpy(SDsum,dev_SDsum,sizeof(float)*nElt1*nstruct/skips,cudaMemcpyDeviceToHost);
      //cout<<i<<" step done"<<endl;
    }

    cudaMemcpy(SD,dev_SD,sizeof(float)*nElt1*SD_store,cudaMemcpyDeviceToHost);
    cudaFree(dev_A);
    cudaFree(dev_B);
    cudaFree(dev_SDsum);
    
    
        
    SDavg=(float *) malloc (sizeof(float)*SD_store);
    for(int j =0; j<SD_store;j++)
    {
      SDavg[j]=0.0;
    }
    cudaMalloc((void **)&dev_SDavg,sizeof(float)*SD_store);
    cudaMemcpy(dev_SDavg,SDavg,sizeof(float)*SD_store,cudaMemcpyHostToDevice);

    dim3 dimBlocks(1024,1,1);
    dim3 dimGrids((SD_store+dimBlocks.x-1)/dimBlocks.x,1,1);
    SDreduce<<<dimGrids,dimBlocks>>>(dev_SD,dev_SDavg,SD_store,nElt1);

    cudaMemcpy(SDavg,dev_SDavg,sizeof(float)*SD_store,cudaMemcpyDeviceToHost);
    cudaFree(dev_SD);
    cudaFree(dev_SDavg);
    
// printing the msd data to be visualized

    FILE *MSDplot=fopen("msd.data","wt");
      fprintf(MSDplot,"# Time (ps) MSD (A^(2)) \n");
    for(int i =1; i<SD_store;i++)
    {
      fprintf(MSDplot,"%f %f\n",(float(i)*timestep)/1000.0,SDavg[i]);
    }
    fclose(MSDplot);
   
    FILE *DIFFplot=fopen("diffco.data","wt");
      fprintf(DIFFplot,"# Time (ps) D (A^(2)/fs) \n");
    for(int i =1; i<SD_store;i++)
    {
      fprintf(DIFFplot,"%f %f\n",(float(i)*timestep)/1000.0,SDavg[i]/(6*i*timestep));
    }
    fclose(DIFFplot);


 
}
