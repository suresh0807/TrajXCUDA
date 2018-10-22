//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"



void compute_VDOS()
{
  
  cout << "metal "<< metal_species<<" counts "<<count_metal<<endl;
  
  int *exch;
  int *dev_exch;
  exch= (int *) malloc (sizeof(int)*nElt1*(nstruct/skip));
 
// Initialize distance matrix and histogram matrix
 for(int l=0; l<nstruct/skip; l++)
    {
      for(int i=0; i< nElt1; i++)
        {
          exch[i+l*nElt1]=0;
        }     
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


    dim3 dimBlock(10,10,10);
    dim3 dimGrid((nElt1+dimBlock.x-1)/dimBlock.x,(count_metal+dimBlock.y-1)/dimBlock.y,((nstruct/skip)+dimBlock.z-1)/dimBlock.z);

//Cuda kernal execution for distance matrix with CUDA timing API commands

    //For surface metal atoms
    cudaMemcpy(dev_exch,exch,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice);
    cudaMemcpy(dev_lattice,lattice,sizeof(int)*6*(nstruct/skip),cudaMemcpyHostToDevice);
    if (cell_type == "orthorhombic") 
    {
    covermat<<<dimGrid,dimBlock>>>(dev_A,dev_B,dev_exch,bondist,nElt1,count_metal,(nstruct/skip),dev_lattice);
    }
    else if (cell_type == "monoclinic")
    {
     covermatmono<<<dimGrid,dimBlock>>>(dev_A,dev_B,dev_exch,bondist,nElt1,count_metal,(nstruct/skip),dev_lattice);
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
  /*  
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

  /*  exchplot.open("exchange.data");
      for(int k=0; k<zsplit; k++)
      {
          exchplot << ztick[k*2+1] <<" "<<exchdensity[k]<<endl;
      }
  
      for(int k=0; k<nElt1; k++)
      {
          exchplot << A[k*3+2] <<" "<<exchsum[k]<<endl;
      }*/
  //  exchplot.close(); 
///////////////////////////////////////////Special atom selection between 2 distances/////////////////////////////////////////////////////////////////////////

  
  if(choose_atoms=="yes")
{
  float surf,avgsurf=0.0;
int dirn;
  if(choose_dirn=="z") dirn=2;
  else if(choose_dirn=="x") dirn=0;
  else if(choose_dirn=="y") dirn=1;
  
cout << Elt1 <<" within "<< minbondist << " and " << maxbondist<<" angstrom of "<< metal_species<<endl;
  
 
  
  int surfatom;
  for(int i=0 ; i<(nstruct/skip); i++)
  {surfatom=0;surf=0.0;
  for(int j=0; j<count_metal;j++)
  {
    if(METAL[j*3+i*count_metal*3+dirn] > mintop && METAL[j*3+i*count_metal*3+dirn] < maxtop)
    {
    surf+=METAL[j*3+i*count_metal*3+dirn];
    surfatom++;
    }
  }
  surf/=surfatom;
  avgsurf+=surf;
  }
  
  avgsurf/=(nstruct/skip);
  
 cout<<"Average surface atoms position in "<<choose_dirn<<" is "<<avgsurf<<endl;
 cout<<"Average number of surface atoms is "<< surfatom<<endl;
  exch= (int *) malloc (sizeof(int)*nElt1*(nstruct/skip));
  Aintnum = (int *) malloc (sizeof(int)*(nstruct/skip));
  Abulknum = (int *) malloc (sizeof(int)*(nstruct/skip));
  
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

    
    cudaMalloc((void **)&dev_A,sizeof(float)*nElt1*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_B,sizeof(float)*count_metal*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_exch,sizeof(int)*nElt1*(nstruct/skip));
    cudaMalloc((void **)&dev_lattice,sizeof(float)*6*(nstruct/skip));
//Copy data from host to device

    cudaMemcpy(dev_A,A,sizeof(float)*nElt1*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_B,METAL,sizeof(float)*count_metal*(nstruct/skip)*3,cudaMemcpyHostToDevice);

//dim3 struct to define elements of the execution configuration


    dim3 dimBlocke(10,10,10);
    dim3 dimGride((nElt1+dimBlocke.x-1)/dimBlocke.x,(count_metal+dimBlocke.y-1)/dimBlocke.y,((nstruct/skip)+dimBlocke.z-1)/dimBlocke.z);

//Cuda kernal execution for distance matrix with CUDA timing API commands

    //For surface metal atoms
    cudaMemcpy(dev_exch,exch,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice);
    cudaMemcpy(dev_lattice,lattice,sizeof(int)*6*(nstruct/skip),cudaMemcpyHostToDevice);
    if (cell_type == "orthorhombic") 
    {
      if(strict_dirn == "yes") {covermat<<<dimGride,dimBlocke>>>(dev_A,dev_exch,minbondist,maxbondist,nElt1,(nstruct/skip),dirn,avgsurf,mintop,maxtop);}
      else {covermat<<<dimGride,dimBlocke>>>(dev_A,dev_B,dev_exch,minbondist,maxbondist,nElt1,count_metal,(nstruct/skip),dev_lattice,dirn,avgsurf,mintop,maxtop);}
    }
    else if (cell_type == "monoclinic")
    {
      if(strict_dirn == "yes") {covermat<<<dimGride,dimBlocke>>>(dev_A,dev_exch,minbondist,maxbondist,nElt1,(nstruct/skip),dirn,avgsurf,mintop,maxtop);} //no lattice needed thus use covermat - no problem
      else {covermatmono<<<dimGride,dimBlocke>>>(dev_A,dev_B,dev_exch,minbondist,maxbondist,nElt1,count_metal,(nstruct/skip),dev_lattice,dirn,avgsurf,mintop,maxtop);}
    }
    
    
    cudaMemcpy(exch,dev_exch,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyDeviceToHost);
    
    
    cudaFree(dev_A);
    cudaFree(dev_B);
    cudaFree(dev_lattice);
    
    //ofstream exchout;
    //exchout.open("exchout.data");
    
    
    for(int l=0; l<nstruct/skip; l++)
      {
	//exchout<<l<<endl;
      for(int i=0; i< nElt1; i++)
        {
	  //exchout<<i<<" "<<exch[i+l*nElt1]<<endl;
          if (exch[i+l*nElt1]>0) {exch[i+l*nElt1]=1; Aintnum[l]++;}
          else if (exch[i+l*nElt1]==0) Abulknum[l]++;
	  //exchout<<i<<" "<<exch[i+l*nElt1]<<endl;
        }
        //exchout<<endl;
      }
      
      //exchout.close();
     cudaMemcpy(dev_exch,exch,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice); 
     
    
}          
   
///////////////////////////////////////////////special selection done - stored in exch///////////////////////////////////////////////////////////////////////////////////

   cout <<"hola"<<endl;
  
    total_time = (nstruct)*timestep;
    cout << "Each frame resolves to "<<timestep<< " fs evolution"<<endl;
    cout << "Total simulation time : " << total_time/1000000 <<" ns"<< endl;
    num_bins =  int ((nstruct)/((diffuse_time*1000)/timestep));
    num_bins -= 1; 
    cout << "The analysis will be done for every "<< diffuse_time <<" ps"<<endl;
    cout <<num_bins <<" number of bins are chosen for the VDOS analysis"<<endl;
    SD_store = int((diffuse_time*1000)/timestep);
    origins = SD_store;
    cout <<"This will correlate for "<<SD_store<<"frames in the input trajectory "<<endl;
    int restarts = origins/skips;
    cout <<restarts<<" restarts from each bin will be taken"<<endl;
    if(Elt1=="all")
    {
    nElt1=natoms;
    nElt2=natoms;
    }
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
    cudaMalloc((void **)&dev_B,sizeof(float)*nstruct*nElt1*3);
    cudaMalloc((void **)&dev_SD,sizeof(float)*nElt1*SD_store*origins/skips);
    cudaMemcpy(dev_A,VEL,sizeof(float)*nstruct*nElt1*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_B,VEL,sizeof(float)*nstruct*nElt1*3,cudaMemcpyHostToDevice);    
    
    cout << "Memory required to store velocity information: "<<nstruct*nElt1*3*2*4/float(1000000000)<< " Gbs " << endl;
    cout << "Memory for local VAF storage: "<<(nElt1*SD_store*restarts)*4/float(1000000000)<< " Gbs" << endl;


    dim3 dimBlocka(32,1,32);
    dim3 dimGrida((nElt1+dimBlocka.x-1)/dimBlocka.x,1,((origins/skips)+dimBlocka.z-1)/dimBlocka.z);

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
      int whichwater;
      if(msd_for=="int"){whichwater=1;}
      else if(msd_for=="bulk"){whichwater=0;}
      else if(msd_for=="all"){whichwater=2;}
      if(choose_atoms=="yes"){whichwater=1;}
cudaMemcpy(dev_SD,SD,sizeof(float)*nElt1*SD_store*origins/skips,cudaMemcpyHostToDevice);
      for(int j=0; j<SD_store; j++)
      {
	//int j=1;
      VAF_calc<<<dimGrida,dimBlocka>>>(dev_A,dev_B,dev_SD,dev_exch,SD_store,nElt1,i,j,origins,skips,whichwater);
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
    
       printf("Padding with zeros and gaussian smoothing... from timteatro\n");
 int padd = 2;
    fftw_complex *VAFavg=(fftw_complex *) fftw_malloc (sizeof(fftw_complex)*SD_store*padd);
    for(int j =0; j<SD_store*padd;j++)
    {
      VAFavg[j]=0.0;
    }    
    for(int j =0; j<SD_store;j++)
    {
      VAFavg[j]=SDavg[j];
    }    
     for(int j =0; j<SD_store;j++)
    {
    //  cout<<j<<" "<<SDavg[j]<<" "<<creal(VAFavg[j])<<endl;
    }
      for(int j =1; j<  SD_store*padd;j++)
    {
      VAFavg[j] /= VAFavg[0];
    }  

    VAFavg[0]=1.0;
    float sigma = float(SD_store) / 2.50;
    
    for(int j =0; j<SD_store*padd;j++)
    {
      VAFavg[j] *=exp( -j * j/(2*sigma*sigma))/(sigma *2.506628274631000);
    } 

     for(int j =1; j<SD_store*padd;j++)
    {
      VAFavg[j] /= VAFavg[0];
    }

    VAFavg[0]=1.0;
// printing the rdf data to be visualized

    fftw_complex    norm = 1.0f;
    fftw_complex   *dft_out;
    fftw_plan       dft_plan;
    
    dft_out = (fftw_complex *) fftw_malloc(SD_store * padd * sizeof(fftw_complex));
    
   dft_plan = fftw_plan_dft_1d(SD_store*padd, VAFavg, dft_out, FFTW_FORWARD, FFTW_ESTIMATE);
   printf("#   done.\n# Executing FFT\n");
   fftw_execute(dft_plan);

   norm = 0.00f;
   for (int m = 0; m <= (SD_store * padd) / 2; m++)
   {
      if (creal(dft_out[m] * conj(dft_out[m])) > creal(norm))
      {
         norm = dft_out[m] * conj(dft_out[m]);
      }
   }
      norm = 1/norm;
   
   
      
    FILE *VAFplot=fopen("vaf.data","wt");
      fprintf(VAFplot,"# Time (ps) VAF (A^(2)/ps^(2)) \n");
    for(int i =0; i<SD_store*padd;i++)
    {

      fprintf(VAFplot,"%f %f\n",(float(i)*timestep)/1000.0,VAFavg[i]);
    }
    fclose(VAFplot);
    
    
    FILE *VDOSplot=fopen("vdos.data","wt");
      fprintf(VDOSplot,"# Time (ps) VAF (A^(2)/ps^(2)) \n");
    for(int i =0; i<(SD_store*padd)/2;i++)
    {

      fprintf(VDOSplot, "%17.9E %17.9E %17.9E %17.9E\n",
              i / (2.99792458E10 * (diffuse_time) * padd  * 1E-12),
              creal(dft_out[i]),
              cimag(dft_out[i]),
              norm * dft_out[i] * conj(dft_out[i]));
    }
    fclose(VDOSplot);    
    
    
}
