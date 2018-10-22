//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"

void compute_angledist(void)
{

  
  //lines for exchange matrix-----------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------------------------
  if(choose_atoms=="yes")
{
cout << Elt1 <<" within "<< bondist <<" angstrom of "<< metal_species<<endl;
  
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
    
    
    for(int l=0; l<nstruct/skip; l++)
      {
      for(int i=0; i< nElt1; i++)
        {
          if (exch[i+l*nElt1]>0) {exch[i+l*nElt1]=1; Aintnum[l]++;}
          else if (exch[i+l*nElt1]==0) Abulknum[l]++;
        }
      }
      
     cudaMemcpy(dev_exch,exch,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice); 
}      
  
//exchange matrix created--------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------
  
  
  
  int *h2onum;
  int *h2oneigh;
  float *h2oxyz;
  
  float *h2odistmat, *dev_h2odistmat;
  
  int OHsamples=4;
  
    h2onum = (int *) malloc (sizeof(int)*nElt1*(nstruct/skip));
    h2oneigh = (int *) malloc (sizeof(int)*OHsamples*nElt1*(nstruct/skip));
    h2oxyz = (float *) malloc (sizeof(float)*13*nElt1*(nstruct/skip));  //Oxyz,H1xyz,H2xyz,Hmidpoint,dipole alignment with surface normal
    h2odistmat = (float *) malloc (sizeof(float)*nElt1*nElt2);
    
  for(int i=0;i<nstruct/skip;i++)
  {
    for(int j=0; j<nElt1; j++)
    {
      h2onum[j+i*nElt1] =0;
      for(int k =0 ;k<OHsamples;k++)
      {
      h2oneigh[j*OHsamples+i*nElt1*OHsamples+k]=0;
      }
      for(int k =0 ;k<9;k++)
      {
      h2oxyz[j*9+i*nElt1*9+k]=0.0;
      }
    }
  }   
      
      
cout <<"gpu begins OH distance computation"<<endl;
//Allocate memory in GPU device

    cudaMalloc((void **)&dev_A,sizeof(float)*nElt1*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_B,sizeof(float)*nElt2*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_h2odistmat,sizeof(float)*nElt1*nElt2);
    cudaMalloc((void **)&dev_lattice,sizeof(float)*6*(nstruct/skip));
//Copy data from host to device

    cudaMemcpy(dev_A,A,sizeof(float)*nElt1*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_B,B,sizeof(float)*nElt2*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_lattice,lattice,sizeof(int)*6*(nstruct/skip),cudaMemcpyHostToDevice);
//dim3 struct to define elements of the execution configuration

    dim3 dimBlock(32,32,1);
    dim3 dimGrid((nElt1+dimBlock.x-1)/dimBlock.x,(nElt2+dimBlock.y-1)/dimBlock.y,1);
    
    bondist=1.10; // water molecule only if OH distance is within this value
    //ofstream neigh;
    //neigh.open("H2O.dat");
    
for(int i = 0; i<(nstruct/skip); i++)
{
  
 // neigh << i<<endl<<endl;
  
   for(int i1=0;i1<nElt1;i1++)
  {
    for(int j1=0; j1<nElt2; j1++)
    {
      h2odistmat[j1+i1*nElt2]=0.0;
    }
  }
  
  cudaMemcpy(dev_h2odistmat,h2odistmat,sizeof(float)*nElt1*nElt2,cudaMemcpyHostToDevice);
//Cuda kernal execution for distance matrix with CUDA timing API commands
    if(cell_type == "orthorhombic")
    {
          Hbondmat<<<dimGrid,dimBlock>>>(dev_A,dev_B,dev_h2odistmat,bondist,nElt1,nElt2,dev_lattice,i);
    }
    else if(cell_type == "monoclinic")
    {
          Hbondmatmono<<<dimGrid,dimBlock>>>(dev_A,dev_B,dev_h2odistmat,bondist,nElt1,nElt2,dev_lattice,i);      
    }

cudaMemcpy(h2odistmat,dev_h2odistmat,sizeof(float)*nElt1*nElt2,cudaMemcpyDeviceToHost);
   
   
   
  for(int i1=0;i1<nElt1;i1++)
  {
    int k=0;
   // neigh << i1<<" ";
    for(int j1=0; j1<nElt2; j1++)
    {
     // neigh <<j1<<" "<<h2odistmat[j1+i1*nElt2]<<" ";
      if(h2odistmat[j1+i1*nElt2] !=0.0)
      {
	h2onum[i1+i*nElt1]++;
	h2oneigh[i1*OHsamples+i*nElt1*OHsamples+k] = j1;
	k++;
      }
    }
    //neigh<<endl;
  }

}   
 
 //neigh.close();
    cudaFree(dev_A);
    cudaFree(dev_B);
    cudaFree(dev_h2odistmat);
    cudaFree(dev_lattice);
    

cout <<"gpu ends"<<endl;
   
cout<<" done !!!"<<endl;    
    
    
    
//Hbondneigh contain the IDs of four possible oxygen neighbour atoms if the hydrogen chosen is participating in a hydrogen bond.
//Onum has the number of oxygen atoms within a specified distance (see input file).

 cout<<"Write Hxyz.dat"<<endl;

  ofstream watercount;
  watercount.open("Hxyz.dat");
  
  
  int *angdist;
  
  angdist = (int *) malloc (sizeof(int)*180);
  
  for(int i=0; i<=180;i++)
  {
    angdist[i]=0;
  }
  
  
for(int i=0; i<(nstruct/skip); i++)
{
 // cout<<i<<endl;
 for(int j=0; j<nElt1;j++)
 {
     
   if(h2onum[j+i*nElt1]==2 && exch[j+i*nElt1]==1)
   //if(h2onum[j+i*nElt1]==2)
   {
        watercount << i <<" "<<j<<" ";
     for(int k=0; k<3;k++)
     {
     h2oxyz[j*9+i*nElt1*9+k]= A[j*3+i*nElt1*3+k];
     watercount << h2oxyz[j*9+i*nElt1*9+k] <<" ";
     }
     for(int k=0; k<3;k++)
     {
     h2oxyz[j*9+i*nElt1*9+3+k]= B[h2oneigh[j*OHsamples+i*nElt1*OHsamples]*3+i*nElt2*3+k];
     //watercount << h2oxyz[j*9+i*nElt1*9+3+k] <<" ";
     }
     for(int k=0; k<3;k++)
     {
     h2oxyz[j*9+i*nElt1*9+6+k]= B[h2oneigh[j*OHsamples+i*nElt1*OHsamples+1]*3+i*nElt2*3+k];
     //watercount << h2oxyz[j*9+i*nElt1*9+6+k] <<" ";
     }
     for(int k=0; k<3;k++)
     {
     h2oxyz[j*9+i*nElt1*9+9+k]=(h2oxyz[j*9+i*nElt1*9+3+k]+h2oxyz[j*9+i*nElt1*9+6+k])/2.0;
     //watercount << h2oxyz[j*9+i*nElt1*9+9+k] <<" ";
     }
     
     h2oxyz[j*9+i*nElt1*9+12]=angle(0,0,0,0,0,1,h2oxyz[j*9+i*nElt1*9],h2oxyz[j*9+i*nElt1*9+1],h2oxyz[j*9+i*nElt1*9+2],h2oxyz[j*9+i*nElt1*9+9],h2oxyz[j*9+i*nElt1*9+9+1],h2oxyz[j*9+i*nElt1*9+9+2]);
     watercount << h2oxyz[j*9+i*nElt1*9+12]*(180.0/3.14159);
     angdist[int(ceil(h2oxyz[j*9+i*nElt1*9+12]*(180.0/3.14159)))]++;
     watercount<<endl;
   }
   
  }
}
  watercount.close();
  

  ofstream angldist;
  angldist.open("Dipole-orient.dat");
  
   for(int i=0; i<=180;i++)
  {
    angldist << i<<" "<<angdist[i]<<endl;
  }
  angldist.close();
cudaFree(dev_exch);
  
}
