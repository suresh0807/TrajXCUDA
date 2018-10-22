//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"

void compute_exchange(void)
{
  
  cout << "metal "<< metal_species<<" counts "<<count_metal<<endl;
  
  int  *exchflip, *dense;
  float *exchsum, *exchdensity;
  exch= (int *) malloc (sizeof(int)*nElt1*(nstruct/skip));
  exchflip= (int *) malloc (sizeof(int)*nElt1*(nstruct/skip));
 
 if(set_max_z == "yes")
{
  minz = set_minz;
  maxz = set_maxz;
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
    

FILE *xyztick = fopen("xyztick.dat","wt");
    for(int i=0;i<zsplit;i++)
    {
      fprintf(xyztick,"%f %f \n",ztick[i*2],ztick[i*2+1]);
    }
fclose(xyztick);

    exchdensity=(float *) malloc (sizeof(float)*zsplit);
    dense = (int *) malloc (sizeof(int)*zsplit);
        for(int k=0;k<zsplit;k++)
        {
          exchdensity[k]=0.0;
	  dense[k]=1;
        }
 
 
// Initialize distance matrix and histogram matrix
 for(int l=0; l<nstruct/skip; l++)
    {
      for(int i=0; i< nElt1; i++)
        {
          exch[i+l*nElt1]=0;
          exchflip[i+l*nElt1]=0;
        }     
      }

 exchsum= (float *) malloc (sizeof(float)*nElt1);   
 
 for(int l=0; l<nElt1; l++)
    {
          exchsum[l]=0.0;
      }
//Allocate memory in GPU device

    cudaMalloc((void **)&dev_B,sizeof(float)*count_metal*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_A,sizeof(float)*nElt1*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_exch,sizeof(int)*nElt1*(nstruct/skip));


//Copy data from host to device

    cudaMemcpy(dev_B,METAL,sizeof(float)*count_metal*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_A,A,sizeof(float)*nElt1*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    

//dim3 struct to define elements of the execution configuration


    dim3 dimBlock(10,10,10);
    dim3 dimGrid((nElt1+dimBlock.x-1)/dimBlock.x,(count_metal+dimBlock.y-1)/dimBlock.y,((nstruct/skip)+dimBlock.z-1)/dimBlock.z);

//Cuda kernal execution for distance matrix with CUDA timing API commands
 
    //For surface metal atoms
    cudaMemcpy(dev_exch,exch,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice);
    if (cell_type == "orthorhombic") 
    {
    covermat<<<dimGrid,dimBlock>>>(dev_A,dev_B,dev_exch,bondist,nElt1,count_metal,(nstruct/skip),xvec,yvec,zvec);
    }
    else if (cell_type == "monoclinic")
    {
     covermatmono<<<dimGrid,dimBlock>>>(dev_A,dev_B,dev_exch,bondist,nElt1,count_metal,(nstruct/skip),xvec,yvec,zvec,yxvec);
    }
    cudaMemcpy(exch,dev_exch,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyDeviceToHost);
    //For substrate atoms
     cudaFree(dev_A);
    cudaFree(dev_B);
    cudaFree(dev_exch);
    

   int startbulk=0, startint=0;
  
    for(int l=0; l<nstruct/skip; l++)
    {
      for(int i=0; i< nElt1; i++)
        {
          if (exch[i+l*nElt1]>0) {exch[i+l*nElt1]=1;}
          exchsum[i]+=exch[i+l*nElt1];
        }     
      }
    
    for(int l=1; l<nstruct/skip; l++)
    {
      for(int i=0; i< nElt1; i++)
        {
          if (exch[i+l*nElt1] != exch[i+l*nElt1-1*nElt1]) {exchflip[i+l*nElt1]=1;} else {exchflip[i+l*nElt1]=0;}
        }     
      }
      int interflip=0, bulkflip=0;
    for(int l=1; l<nstruct/skip; l++)
    {
      for(int i=0; i< nElt1; i++)
        {
          if (exch[i]==1 && exchflip[i+l*nElt1]==1) {interflip++;} 
          else if (exch[i]==0 && exchflip[i+l*nElt1]==1){bulkflip++;}
        }     
      } 
      for(int i=0; i< nElt1; i++)
        {
          if (exch[i]==1) {startint++;}
          else if (exch[i]==0){startbulk++;}
        }
      cout <<interflip/startint<<"  "<<bulkflip/startbulk<<endl;
    
    cout<<"Its okay"<<endl;
    
        for(int k=0;k<zsplit;k++)
        {
	  for(int i=0; i< nElt1; i++)
	  {
	    if(A[i*3+2] > ztick[k*2] && A[i*3+2] < ztick[k*2+1])
	    {
	     dense[k]++;
             exchdensity[k]+=exchsum[i];
	    }
	  }
	  exchdensity[k]=exchdensity[k]/float(dense[k]);
        } 

      
      
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
  
    /*  for(int k=0; k<nElt1; k++)
      {
          exchplot << A[k*3+2] <<" "<<exchsum[k]<<endl;
      }*/
  //  exchplot.close(); 
  free(exchsum);
  free(exch);
  free(A);
  free(METAL);
  
}
