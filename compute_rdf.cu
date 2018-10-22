//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"
float *rdf;
float *dev_rdf;
int bin;
float bin_rad;
float pi;
float *rdfavgf;
float rho;
float *dev_rdfavgf;
float *rdfavg;
float *dev_rdfavg;


void compute_rdf(void)
{
 /*
  if(multiply=="yes")
  {
  float *Aext, *Bext;
  
  nElt1ext= nElt1*4;
  nElt2ext= nElt2*4;
  
  Aext= (float *) malloc (sizeof(float)*nElt1*3*(nstruct/skip));
  Bext= (float *) malloc (sizeof(float)*nElt2*3*(nstruct/skip));
  
  for(int i=0; i<(nstruct/skip); i++)
  {
    for(int j=0;j<nElt1;j++)
    {
      Aext[j*3+ i*nElt1ext*3]=A[j*3+ i*nElt1*3]
    }
    for(int j=nElt1;j<nElt1ext/3;j++)
    {
      Aext[j*3+ i*nElt1ext*3]=A[(j-nElt1)*3+ i*nElt1*3]+lattice[(j-nElt1)*6+i*nElt1*6];
    }
    for(int j=nElt1;j<nElt1ext/2;j++)
    {
      Aext[j*3+ i*nElt1ext*3]  =A[(j-(nElt1*2))*3+ i*nElt1*3]  +lattice[(j-(nElt1*2))*6+i*nElt1*6+3];
      Aext[j*3+ i*nElt1ext*3+1]=A[(j-(nElt1*2))*3+ i*nElt1*3+1]+lattice[(j-(nElt1*2))*6+i*nElt1*6+3+1];
    }
    for(int j=nElt1*2;j<nElt1ext;j++)
    {
      Aext[j*3+ i*nElt1ext*3]  =A[(j-(nElt1*3))*3+ i*nElt1*3]  +lattice[(j-(nElt1*3))*6+i*nElt1*6]  +lattice[(j-(nElt1*3))*6+i*nElt1*6+3];
      Aext[j*3+ i*nElt1ext*3+1]=A[(j-(nElt1*3))*3+ i*nElt1*3+1]+lattice[(j-(nElt1*3))*6+i*nElt1*6+1];
    }
  }
  
  }
  */
  int whichwater;
  int *exch1, *exch2;
  int *dev_exch1, *dev_exch2;
  float *dev_A1, *dev_A2;
  
  exch1= (int *) malloc (sizeof(int)*nElt1*(nstruct/skip));
  exch2= (int *) malloc (sizeof(int)*nElt2*(nstruct/skip));

  Aintnum = (int *) malloc (sizeof(int)*(nstruct/skip));
  Bintnum = (int *) malloc (sizeof(int)*(nstruct/skip));
  Abulknum = (int *) malloc (sizeof(int)*(nstruct/skip));
  Bbulknum = (int *) malloc (sizeof(int)*(nstruct/skip));
// Initialize distance matrix and histogram matrix
 for(int l=0; l<nstruct/skip; l++)
    {
      for(int i=0; i< nElt1; i++)
        {
          exch1[i+l*nElt1]=0;
        }
      for(int i=0; i< nElt2; i++)
        {
          exch2[i+l*nElt2]=0;
        }
        Abulknum[l]=0;
        Bbulknum[l]=0;
        Aintnum[l]=0;
        Bintnum[l]=0;
      }
  
      
      float GLOBALAminz=A[dircover];
      float GLOBALAmaxz=A[dircover];
      float MIDAz=0.0, MINAz=0.0, MAXAz=0.0;
      
  for(int i=0; i<1; i++)
  {
    for(int j=1; j<nElt1; j++)
    {
    if(A[j*3+i*nElt1*3+dircover] > GLOBALAmaxz)    GLOBALAmaxz=A[j*3+i*nElt1*3+dircover];
    else if(A[j*3+i*nElt1*3+dircover] < GLOBALAminz)    GLOBALAminz=A[j*3+i*nElt1*3+dircover];
    }
  }
  
  MIDAz=(GLOBALAmaxz+GLOBALAminz)/2.0;
  
  MAXAz=MIDAz+(bondist/2.0);
  MINAz=MIDAz-(bondist/2.0);
     
  if(rdf_between == "yes") {MAXAz=maxtop; MINAz=mintop; MIDAz=(MAXAz+MINAz)/2.0;bondist=MAXAz-MINAz;}
  
  cudaMalloc((void **)&dev_exch1,sizeof(int)*nElt1*(nstruct/skip));
  cudaMalloc((void **)&dev_exch2,sizeof(int)*nElt2*(nstruct/skip));
  cudaMemcpy(dev_exch1,exch1,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice);
  cudaMemcpy(dev_exch2,exch2,sizeof(int)*nElt2*(nstruct/skip),cudaMemcpyHostToDevice);
  cudaMalloc((void **)&dev_lattice,sizeof(float)*6*(nstruct/skip));
  cudaMemcpy(dev_lattice,lattice,sizeof(int)*6*(nstruct/skip),cudaMemcpyHostToDevice);
  
  
  if(rdf_for=="int")
  {
    
  cout << "metal "<< metal_species<<" counts "<<count_metal<<endl;

//Allocate memory in GPU device

    cudaMalloc((void **)&dev_B,sizeof(float)*count_metal*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_A1,sizeof(float)*nElt1*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_A2,sizeof(float)*nElt2*(nstruct/skip)*3);
    

//Copy data from host to device

    cudaMemcpy(dev_B,METAL,sizeof(float)*count_metal*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_A1,A,sizeof(float)*nElt1*(nstruct/skip)*3,cudaMemcpyHostToDevice);


//dim3 struct to define elements of the execution configuration

    dim3 dimBlock(10,10,10);
    dim3 dimGrid((nElt1+dimBlock.x-1)/dimBlock.x,(count_metal+dimBlock.y-1)/dimBlock.y,((nstruct/skip)+dimBlock.z-1)/dimBlock.z);

//Cuda kernal execution for distance matrix with CUDA timing API commands

  
    
    
    if (cell_type == "orthorhombic")
    {
    covermat<<<dimGrid,dimBlock>>>(dev_A1,dev_B,dev_exch1,bondist,nElt1,count_metal,(nstruct/skip),dev_lattice);
    }
    else if (cell_type == "monoclinic")
    {
     covermatmono<<<dimGrid,dimBlock>>>(dev_A1,dev_B,dev_exch1,bondist,nElt1,count_metal,(nstruct/skip),dev_lattice);
    }
    cudaMemcpy(exch1,dev_exch1,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyDeviceToHost);
    //For substrate atoms
     cudaFree(dev_A1);

    dim3 dimBlock2(10,10,10);
    dim3 dimGrid2((nElt2+dimBlock2.x-1)/dimBlock2.x,(count_metal+dimBlock2.y-1)/dimBlock2.y,((nstruct/skip)+dimBlock2.z-1)/dimBlock2.z);

//Cuda kernal execution for distance matrix with CUDA timing API commands

    cudaMemcpy(dev_A2,B,sizeof(float)*nElt2*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    
    if (cell_type == "orthorhombic")
    {
    covermat<<<dimGrid2,dimBlock2>>>(dev_A2,dev_B,dev_exch2,bondist,nElt2,count_metal,(nstruct/skip),dev_lattice);
    }
    else if (cell_type == "monoclinic")
    {
     covermatmono<<<dimGrid2,dimBlock2>>>(dev_A2,dev_B,dev_exch2,bondist,nElt2,count_metal,(nstruct/skip),dev_lattice);
    }
    cudaMemcpy(exch2,dev_exch2,sizeof(int)*nElt2*(nstruct/skip),cudaMemcpyDeviceToHost);
    //For substrate atoms
     cudaFree(dev_A2);
     cudaFree(dev_B);
    //cudaFree(dev_exch);
     
     
    for(int l=0; l<nstruct/skip; l++)
    {
      for(int i=0; i< nElt1; i++)
        {
          if (exch1[i+l*nElt1]>0) {exch1[i+l*nElt1]=1; Aintnum[l]++;}
          else if (exch1[i+l*nElt1]==0) Abulknum[l]++;
        }
      for(int i=0; i< nElt2; i++)
        {
          if (exch2[i+l*nElt2]>0) {exch2[i+l*nElt2]=1; Bintnum[l]++;}
          else if (exch2[i+l*nElt2]==0) Bbulknum[l]++;
        }
      }
    cudaMemcpy(dev_exch1,exch1,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice);
    cudaMemcpy(dev_exch2,exch2,sizeof(int)*nElt2*(nstruct/skip),cudaMemcpyHostToDevice);

    
      whichwater=1;
  }
  
    
  else if(rdf_for=="bulk")
  {
    
//Allocate memory in GPU device


    cudaMalloc((void **)&dev_A1,sizeof(float)*nElt1*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_A2,sizeof(float)*nElt2*(nstruct/skip)*3);

    
//Copy data from host to device

    cudaMemcpy(dev_A1,A,sizeof(float)*nElt1*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    
    //dim3 struct to define elements of the execution configuration

    dim3 dimBlock(32,1,32);
    dim3 dimGrid((nElt1+dimBlock.x-1)/dimBlock.x,1,((nstruct/skip)+dimBlock.z-1)/dimBlock.z);

//Cuda kernal execution for distance matrix with CUDA timing API commands

  
    if (cell_type == "orthorhombic")
    {
    covermat<<<dimGrid,dimBlock>>>(dev_A1,dev_exch1,nElt1,(nstruct/skip),MINAz,MAXAz,dircover);
    }
    else if (cell_type == "monoclinic")
    {
     covermatmono<<<dimGrid,dimBlock>>>(dev_A1,dev_exch1,nElt1,(nstruct/skip),MINAz,MAXAz,dircover);
    }
    cudaMemcpy(exch1,dev_exch1,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyDeviceToHost);
    //For substrate atoms
     cudaFree(dev_A1);

     cudaMemcpy(dev_A2,B,sizeof(float)*nElt2*(nstruct/skip)*3,cudaMemcpyHostToDevice);
     
    dim3 dimBlock2(32,1,32);
    dim3 dimGrid2((nElt2+dimBlock2.x-1)/dimBlock2.x,1,((nstruct/skip)+dimBlock2.z-1)/dimBlock2.z);

//Cuda kernal execution for distance matrix with CUDA timing API commands

    if (cell_type == "orthorhombic")
    {
    covermat<<<dimGrid2,dimBlock2>>>(dev_A2,dev_exch2,nElt2,(nstruct/skip),MINAz,MAXAz,dircover);
    }
    else if (cell_type == "monoclinic")
    {
     covermatmono<<<dimGrid2,dimBlock2>>>(dev_A2,dev_exch2,nElt2,(nstruct/skip),MINAz,MAXAz,dircover);
    }
    cudaMemcpy(exch2,dev_exch2,sizeof(int)*nElt2*(nstruct/skip),cudaMemcpyDeviceToHost);
    //For substrate atoms
     cudaFree(dev_A2);
    //cudaFree(dev_exch);
     
    
    for(int l=0; l<nstruct/skip; l++)
    {
      for(int i=0; i< nElt1; i++)
        {
          if (exch1[i+l*nElt1]>0) {exch1[i+l*nElt1]=1; Abulknum[l]++;}
          else if (exch1[i+l*nElt1]==0) Aintnum[l]++;
        }
      for(int i=0; i< nElt2; i++)
        {
          if (exch2[i+l*nElt2]>0) {exch2[i+l*nElt2]=1; Bbulknum[l]++;}
          else if (exch2[i+l*nElt2]==0) Bintnum[l]++;
        }
      }
    cudaMemcpy(dev_exch1,exch1,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice);
    cudaMemcpy(dev_exch2,exch2,sizeof(int)*nElt2*(nstruct/skip),cudaMemcpyHostToDevice);

      
      whichwater=1;
    
  }

  
  if(rdf_for=="all"){whichwater=2;}
  
pi=3.14159;
    rdf= (float *) malloc (sizeof(float)*nElt1*bin*(nstruct/skip));
// Initialize distance matrix and histogram matrix
 for(int l=0; l<nstruct/skip; l++)
    {
      for(int i=0; i< nElt1; i++)
      {
        for(int j=0; j< bin; j++)
        {
          rdf[j+i*bin+l*nElt1*bin]=0.0;
        }
      }
    }

//Allocate memory in GPU device

    cudaMalloc((void **)&dev_A,sizeof(float)*nElt1*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_B,sizeof(float)*nElt2*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_rdf,sizeof(float)*nElt1*bin*(nstruct/skip));
    double evvalo = double(((nElt1+nElt2)/1000.0*(nstruct/skip)/1000.0*3)*4/1000.0);
    cout << evvalo  \
         << "Gbs of memory needed to store coordinates" <<endl;
             evvalo = double(nElt1/1000.0)*double(nstruct/skip)*double(bin/1000.0)*(4.0/1000.0);
    std::cout << evvalo  \
         << "Gbs of memory needed to store distance matrix" <<endl;	 
    
//Copy data from host to device

    cudaMemcpy(dev_A,A,sizeof(float)*nElt1*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_B,B,sizeof(float)*nElt2*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_rdf,rdf,sizeof(float)*nElt1*bin*(nstruct/skip),cudaMemcpyHostToDevice);

//dim3 struct to define elements of the execution configuration
      


    //bin_rad=(xvec/2)/float(bin);
    bin_rad = rdf_max_rad/float(bin);
    dim3 dimBlockrdf(10,10,10);
    dim3 dimGridrdf((nElt1+dimBlockrdf.x-1)/dimBlockrdf.x,(nElt2+dimBlockrdf.y-1)/dimBlockrdf.y,((nstruct/skip)+dimBlockrdf.z-1)/dimBlockrdf.z);

//Cuda kernal execution for distance matrix with CUDA timing API commands
    
    
    
 //rdf kernels are different for int, bulk and all. For all, usual rdf is done and it makes no sense to compute this for interfacial systems.
 //For bulk, a mid point is defined and a box is cut out from the middle and for all atoms
 //in the box, every other atom within the bondist is correlarted.
 //For int, a small region from the surface, effectively the first layer is chosen.
 //and a hemisphere is sampled for each interface, giving the final rdf.
    if(cell_type == "orthorhombic")
    {
      if(rdf_for =="int")
      {
	distmat<<<dimGridrdf,dimBlockrdf>>>(dev_A,dev_B,dev_rdf,dev_exch1,dev_exch2,bin_rad,bin,nElt1,nElt2,(nstruct/skip),dev_lattice,rdf_max_rad,whichwater,MIDAz,dircover);
      }
      else
    distmat<<<dimGridrdf,dimBlockrdf>>>(dev_A,dev_B,dev_rdf,dev_exch1,dev_exch2,bin_rad,bin,nElt1,nElt2,(nstruct/skip),dev_lattice,rdf_max_rad,whichwater);
    }
    else if(cell_type == "monoclinic")
    {
      if(rdf_for =="int")
      {
	distmatmono<<<dimGridrdf,dimBlockrdf>>>(dev_A,dev_B,dev_rdf,dev_exch1,dev_exch2,bin_rad,bin,nElt1,nElt2,(nstruct/skip),dev_lattice,rdf_max_rad,whichwater,MIDAz,dircover);
      }
      else
    distmatmono<<<dimGridrdf,dimBlockrdf>>>(dev_A,dev_B,dev_rdf,dev_exch1,dev_exch2,bin_rad,bin,nElt1,nElt2,(nstruct/skip),dev_lattice,rdf_max_rad,whichwater);
    }
    
    cudaMemcpy(rdf,dev_rdf,sizeof(float)*nElt1*bin*(nstruct/skip),cudaMemcpyDeviceToHost);

    cudaFree(dev_A);
    cudaFree(dev_B);
    cudaFree(dev_exch1);
    cudaFree(dev_exch2);
    cudaFree(dev_lattice);
    
    
    rdfavgf=(float *) malloc (sizeof(float)*bin*(nstruct/skip));
    for(int i=0; i<(nstruct/skip); i++)
    {
      for(int j =0; j<bin;j++)
      {
        rdfavgf[j+i*bin]=0.0;
      }
    }
    cudaMalloc((void **)&dev_rdfavgf,sizeof(float)*bin*(nstruct/skip));
    cudaMemcpy(dev_rdfavgf,rdfavgf,sizeof(float)*bin*(nstruct/skip),cudaMemcpyHostToDevice);

    dim3 dimBlockf(32,32,1);
    dim3 dimGridf((bin+dimBlockf.x-1)/dimBlockf.x,((nstruct/skip)+dimBlockf.y-1)/dimBlockf.y,1);
    reducef<<<dimGridf,dimBlockf>>>(dev_rdf,dev_rdfavgf,bin,(nstruct/skip),nElt1);
    cudaMemcpy(rdfavgf,dev_rdfavgf,sizeof(float)*bin*(nstruct/skip),cudaMemcpyDeviceToHost);
    cudaFree(dev_rdf);

    rdfavg=(float *) malloc (sizeof(float)*bin);
    for(int j =0; j<bin;j++)
    {
      rdfavg[j]=0.0;
    }
    cudaMalloc((void **)&dev_rdfavg,sizeof(float)*bin);
    cudaMemcpy(dev_rdfavg,rdfavg,sizeof(float)*bin,cudaMemcpyHostToDevice);

    dim3 dimBlocks(1024,1,1);
    dim3 dimGrids((bin+dimBlocks.x-1)/dimBlocks.x,1,1);
    reduce<<<dimGrids,dimBlocks>>>(dev_rdfavgf,dev_rdfavg,bin,(nstruct/skip));

    cudaMemcpy(rdfavg,dev_rdfavg,sizeof(float)*bin,cudaMemcpyDeviceToHost);
    cudaFree(dev_rdfavgf);
    cudaFree(dev_rdfavg);

    for(int k=0; k< bin; k++)
        {
          cout<<rdfavg[k]<<" ";
        }
        cout<<endl;
    
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
    
    float AVGAintnum=0.0, AVGAbulknum=0.0;
    float AVGBintnum=0.0, AVGBbulknum=0.0;
    for(int l=0; l<nstruct/skip; l++)
    {
        AVGAintnum+=Aintnum[l];
        AVGAbulknum+=Abulknum[l];
	AVGBintnum+=Bintnum[l];
        AVGBbulknum+=Bbulknum[l];
    }
AVGAintnum/=float(nstruct/skip);
AVGAbulknum/=float(nstruct/skip);
AVGBintnum/=float(nstruct/skip);
AVGBbulknum/=float(nstruct/skip);
cout <<AVGAintnum<<" "<<AVGAbulknum<<" "<<AVGBintnum<<" "<<AVGBbulknum<<endl;
// printing the rdf data to be visualized
    rho = 0.0;
    
    float cellvolume,intvolume,bulkvolume;
    
    if(cell_type=="orthorhombic")
    {
      cellvolume = xvec*yvec*zvec;
    if(vacuum.at(0) =='z') {intvolume = xvec*yvec*(MIDAz-GLOBALAminz);bulkvolume = xvec*yvec*bondist;}
    else if(vacuum.at(0) =='x') {intvolume = zvec*yvec*(MIDAz-GLOBALAminz);bulkvolume = zvec*yvec*bondist;}
    else if(vacuum.at(0) =='y') {intvolume = xvec*zvec*(MIDAz-GLOBALAminz);bulkvolume = xvec*zvec*bondist;}
    //intvolume = xvec*yvec*(MIDAz-GLOBALAminz);
    //intvolume = xvec*yvec*bondist; // effectively considered as only one interface made of two hemispheres from each interfacial region 
    //bulkvolume = xvec*yvec*bondist;
    }
    else if(cell_type=="monoclinic")
    {
    cellvolume= xvec*xvec*zvec*0.866025403; //*sin(60) for rhombohedral 
    intvolume = xvec*xvec*(MIDAz-GLOBALAminz)*0.866025403;
    //intvolume=xvec*xvec*bondist*0.866025403;
    bulkvolume=xvec*xvec*bondist*0.866025403;
    }
    
    
    if (Elt1 != rdf_metal_exclude && Elt2 != rdf_metal_exclude && rdf_for == "all")
    {
      rho = (nElt1*nElt2)/(cellvolume - (float(count_metal) * metal_atom_volume));
      cout<<(nElt1*nElt2)<<" "<<cellvolume<<" "<<(float(count_metal) * metal_atom_volume)<<" "<<rho<<endl;
    }
    else if (Elt1 != rdf_metal_exclude && Elt2 != rdf_metal_exclude && rdf_for == "int")
    {
      for(int l=0;l<(nstruct/skip);l++)
      {
      //rho += (Aintnum[l]*Bintnum[l])/(intvolume); // since, two interfacial regions
      rho += (Aintnum[l]*nElt2/2.0)/(intvolume);
      }
      rho /=float(nstruct/skip);
      rho/=2.0; // since only the surface area of hemisphere taken for normalizing
      cout<<(AVGAintnum*AVGBintnum)/intvolume<<" "<<rho<<endl;
    }
    else if (Elt1 != rdf_metal_exclude && Elt2 != rdf_metal_exclude && rdf_for == "bulk")
    {
      for(int l=0;l<(nstruct/skip);l++)
      {  
      rho += (Abulknum[l]*Bbulknum[l])/bulkvolume;
      }
      rho /=float(nstruct/skip);
      cout<<(AVGAbulknum*AVGBbulknum)<<" "<< bulkvolume <<" "<<(AVGAbulknum*AVGBbulknum)/bulkvolume<<" "<<rho<<endl;
    }
    else  rho=(nElt1*nElt2)/cellvolume ;
    
    FILE *rdfplot=fopen("rdf.data","wt");
    for(int i =1; i<bin;i++)
    {
      fprintf(rdfplot,"%f %f\n",((bin_rad*(i+1))+(bin_rad*(i)))/2,((rdfavg[i])/(nstruct/skip))/(rho*4*pi*(square(bin_rad*(i+1)))*bin_rad));
      //fprintf(rdfplot,"%f %f\n",((bin_rad*(i+1))+(bin_rad*(i)))/2,((rdfavg[i])/(nstruct/skip)));
      //cout<<i<<" "<<rho<<" "<<4*pi*(square(bin_rad*(i+1)))*bin_rad<<" "<<rho*4*pi*(square(bin_rad*(i+1)))*bin_rad<<endl;
    }
    fclose(rdfplot);
  
}
