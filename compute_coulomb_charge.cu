//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"

void compute_coulomb_charge(void)
{


  float *coulomb_distmat;
  float *coulomb_final;
  coulomb_distmat=(float *) malloc (sizeof(float)*nElt1*nElt2*(nstruct/skip));
  coulomb_final=(float *) malloc (sizeof(float)*(nstruct/skip));
  for(int i=0;i<(nstruct/skip);i++)
  {
    coulomb_final[i]=0.0;
  for(int j=0;j<nElt1;j++)
  {
    
  for(int k=0;k<nElt2;k++)
  {

    coulomb_distmat[k+j*nElt2+i*nElt1*nElt2]=0.0;

  }
  }
  }
  float *dev_coulomb_distmat;
  
  
  
  
  cudaMalloc((void **)&dev_lattice,sizeof(float)*6*(nstruct/skip));
  cudaMalloc((void **)&dev_A,sizeof(float)*nElt1*(nstruct/skip)*4);
  cudaMalloc((void **)&dev_B,sizeof(float)*nElt2*(nstruct/skip)*4);
  cudaMalloc((void **)&dev_coulomb_distmat,sizeof(float)*nElt1*nElt2*(nstruct/skip));
    
//Copy data from host to device

  cudaMemcpy(dev_A,A_all,sizeof(float)*nElt1*(nstruct/skip)*4,cudaMemcpyHostToDevice);
  cudaMemcpy(dev_B,B_all,sizeof(float)*nElt2*(nstruct/skip)*4,cudaMemcpyHostToDevice);
  cudaMemcpy(dev_lattice,lattice,sizeof(float)*6*(nstruct/skip),cudaMemcpyHostToDevice);  
  cudaMemcpy(dev_coulomb_distmat,coulomb_distmat,sizeof(float)*nElt1*nElt2*(nstruct/skip),cudaMemcpyHostToDevice);
  
  cout<<"GPU coulomb matrix computation"<<endl;
  
    dim3 dimBlock(10,10,10);
    dim3 dimGrid((nElt1+dimBlock.x-1)/dimBlock.x,(nElt2+dimBlock.y-1)/dimBlock.y,((nstruct/skip)+dimBlock.z-1)/dimBlock.z);
   
    
    if(cell_type == "orthorhombic")
    {
	distmat<<<dimGrid,dimBlock>>>(dev_A,dev_B,dev_coulomb_distmat,nElt1,nElt2,(nstruct/skip),dev_lattice,Coulomb_cutoff);
    }
    else if(cell_type == "monoclinic")
    {
	distmatmono<<<dimGrid,dimBlock>>>(dev_A,dev_B,dev_coulomb_distmat,nElt1,nElt2,(nstruct/skip),dev_lattice,Coulomb_cutoff);
    }
    
    cudaMemcpy(coulomb_distmat,dev_coulomb_distmat,sizeof(float)*nElt1*nElt2*(nstruct/skip),cudaMemcpyDeviceToHost);
    
    
    cout<<"GPU computation completed: freeing memory"<<endl;
    
    cudaFree(dev_A);
    cudaFree(dev_B);
    cudaFree(dev_coulomb_distmat);
    cudaFree(dev_lattice);
    
    for(int i=0;i<(nstruct/skip);i++)
    {
    for(int j=0;j<nElt1;j++)
    {
    for(int k=0;k<nElt2;k++)
    {
    coulomb_final[i]+=coulomb_distmat[k+j*nElt2+i*nElt1*nElt2];
    }
    }
    }
  
  

    ofstream coulomb;
    coulomb.open("coulomb_E.dat");
    for(int i=0;i<(nstruct/skip);i++)
    {
    coulomb<<i<<" "<<coulomb_final[i]/2<<endl;
    }
    coulomb.close();
    
    //free(A_all);
    //free(B_all);
  }





