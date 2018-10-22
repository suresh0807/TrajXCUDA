//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"

void choose_atoms()
{
  
   
  cout << Elt1 <"within "<< bondist <<" \AA of"<< metal_species<<endl;
  
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
    }
    cudaMemcpy(exch,dev_exch,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyDeviceToHost);
    
    
    cudaFree(dev_A);
    cudaFree(dev_B);
    cudaFree(dev_lattice);
    cudaFree(dev_exch);
    
    for(int l=0; l<nstruct/skip; l++)
      {
      for(int i=0; i< nElt1; i++)
        {
          if (exch[i+l*nElt1]>0) {exch[i+l*nElt1]=1; Aintnum[l]++;}
          else if (exch[i+l*nElt1]==0) Abulknum[l]++;
        }
      }

          
}