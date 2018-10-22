//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"




void compute_diffusion()
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
//cout<<"hi "<< (nstruct/skip)*nElt1*3<<endl;
      cudaMemcpy(dev_A,A,sizeof(float)*(nstruct/skip)*nElt1*3,cudaMemcpyHostToDevice);
//cout<<"hi"<<endl;
    cudaMemcpy(dev_B,B,sizeof(float)*(nstruct/skip)*nElt2*3,cudaMemcpyHostToDevice);
        cout << "Memory required to store coordinate information: "<<nstruct*nElt1*3*2*4/float(1000000000)<< " Gbs " << endl;
    cout << "Memory for local MSD storage: "<<((nElt1*nstruct/skips)+(nElt1*SD_store))*4/float(1000000000)<< " Gbs" << endl;
    }

    
   // cout << "Memory required to store coordinate information: "<<nstruct*nElt1*3*2*4/float(1000000000)<< " Gbs " << endl;
   // cout << "Memory for local MSD storage: "<<nElt1*SD_store*origins/skips*4/float(1000000000)<< " Gbs" << endl;


    dim3 dimBlock(32,1,32);
    dim3 dimGrid((nElt1+dimBlock.x-1)/dimBlock.x,1,((origins/skips)+dimBlock.z-1)/dimBlock.z);

    //dim3 dimBlock(10,10,10);
    //dim3 dimGrid((nElt1+dimBlock.x-1)/dimBlock.x,(SD_store+dimBlock.y-1)/dimBlock.y,((origins/skips)+dimBlock.z-1)/dimBlock.z);    
    
  //    cudaMemcpy(dev_SD,SD,sizeof(float)*nElt1*SD_store*origins,cudaMemcpyHostToDevice);    

    for(int i=0; i<num_bins;i++)
    //for(int i=0; i<num_bins;i++)
    {
             for(int i1=0;i1<origins/skips;i1++)
              {
               for(int j1=0;j1<SD_store;j1++)
                {
	         for(int k1=0; k1<nElt1; k1++)
	          {
                    SD[k1+j1*nElt1+i1*nElt1*SD_store]=0.0;
	          }
                }
              }
      
cudaMemcpy(dev_SD,SD,sizeof(float)*nElt1*SD_store*origins/skips,cudaMemcpyHostToDevice);
      for(int j=1; j<SD_store; j++)
      {
	//int j=1;
      SD_calc<<<dimGrid,dimBlock>>>(dev_A,dev_B,dev_SD,SD_store,nElt1,i,j,origins,skips,xsrt,xend,xski);
      }
cudaMemcpy(SD,dev_SD,sizeof(float)*nElt1*SD_store*origins/skips,cudaMemcpyDeviceToHost);
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
    cudaFree(dev_B);
    cudaFree(dev_SD);
     
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

    
    SDavg=(float *) malloc (sizeof(float)*SD_store);
    for(int j =0; j<SD_store;j++)
    {
      SDavg[j]=0.0;
    }
    cudaMalloc((void **)&dev_SDavg,sizeof(float)*SD_store);
    cudaMemcpy(dev_SDavg,SDavg,sizeof(float)*SD_store,cudaMemcpyHostToDevice);

    dim3 dimBlocks(1024,1,1);
    dim3 dimGrids((SD_store+dimBlocks.x-1)/dimBlocks.x,1,1);
    SDreduce<<<dimGrids,dimBlocks>>>(dev_SDavgf,dev_SDavg,SD_store,nElt1);

    cudaMemcpy(SDavg,dev_SDavg,sizeof(float)*SD_store,cudaMemcpyDeviceToHost);
    
    cudaFree(dev_SDavgf);
    cudaFree(dev_SDavg);
    
// printing the rdf data to be visualized

    FILE *MSDplot=fopen("msd.data","wt");
      fprintf(MSDplot,"# Time (ps) MSD (A^(2)) \n");
    for(int i =1; i<SD_store;i++)
    {
      fprintf(MSDplot,"%f %f\n",(float(i)*timestep)/1000.0,SDavg[i]/num_bins);
    }
    fclose(MSDplot);
   
    FILE *DIFFplot=fopen("diffco.data","wt");
      fprintf(DIFFplot,"# Time (ps) D (A^(2)/fs) \n");
    for(int i =1; i<SD_store;i++)
    {
      fprintf(DIFFplot,"%f %f\n",(float(i)*timestep)/1000.0,SDavg[i]/(num_bins*6*i*timestep));
    }
    fclose(DIFFplot);


 
}
