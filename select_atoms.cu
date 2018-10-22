//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"



void select_atoms(string whichmat)
{
  
  if(whichmat != "box")
  {
    
    
    cout << "metal "<< metal_species<<" counts "<<count_metal<<endl;
  
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
  
  
  if(msd_for=="int")
  {
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
          if (exch[i+l*nElt1]>0) {exch[i+l*nElt1]=1; Aintnum[l]++;}
          else if (exch[i+l*nElt1]==0) Abulknum[l]++;
        }
      }
    cudaMemcpy(dev_exch,exch,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice);
    
      ofstream exchplot;
      }
  
   
  
  else if(msd_for=="bulk")
  {
    //Allocate memory in GPU device

    cudaMalloc((void **)&dev_A,sizeof(float)*nElt1*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_exch,sizeof(int)*nElt1*(nstruct/skip));
    
//Copy data from host to device

    
    cudaMemcpy(dev_A,A,sizeof(float)*nElt1*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    

//dim3 struct to define elements of the execution configuration


    dim3 dimBlocke(32,1,32);
    dim3 dimGride((nElt1+dimBlocke.x-1)/dimBlocke.x,1,((nstruct/skip)+dimBlocke.z-1)/dimBlocke.z);

//Cuda kernal execution for distance matrix with CUDA timing API commands

    //For surface metal atoms
    cudaMemcpy(dev_exch,exch,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice);
    if (cell_type == "orthorhombic") 
    {
    covermat<<<dimGride,dimBlocke>>>(dev_A,dev_exch,nElt1,(nstruct/skip),MINAz,MAXAz,dircover);
    }
    else if (cell_type == "monoclinic")
    {
    covermatmono<<<dimGride,dimBlocke>>>(dev_A,dev_exch,nElt1,(nstruct/skip),MINAz,MAXAz,dircover);
      //covermat<<<dimGride,dimBlocke>>>(dev_A,dev_B,dev_exch,bondist,nElt1,count_metal,(nstruct/skip),dev_lattice);
    }
    cudaMemcpy(exch,dev_exch,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyDeviceToHost);
    //For substrate atoms
     cudaFree(dev_A);
    //cudaFree(dev_exch);
      for(int l=0; l<nstruct/skip; l++)
      {
      for(int i=0; i< nElt1; i++)
        {
          if (exch[i+l*nElt1]>0) {exch[i+l*nElt1]=1; Abulknum[l]++;}
          else if (exch[i+l*nElt1]==0) Aintnum[l]++;
        }
      }
    cudaMemcpy(dev_exch,exch,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice);
    ofstream exchplot;
 

  }
  
  else if (msd_for=="all")
  {
    cudaMalloc((void **)&dev_exch,sizeof(int)*nElt1*(nstruct/skip));
    cudaMemcpy(dev_exch,exch,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice);
    
    MSDchaos=(int *) malloc (sizeof(int)*nElt1*(nstruct/skip));
    for(int l=0; l<(nstruct/skip);l++)
    {
    for(int i=0;i<nElt1;i++)
    {
      MSDchaos[i+l*nElt1]=0;
    }
    }
    
    //MSDsteps=ceil((GLOBALAmaxz-GLOBALAminz))/2;
    //MSDsteps=10;
    cout<<"The cell will be divided into "<<MSDsteps<<" slices along vacuum direction for MSD computation"<<endl;
    MSDint = (GLOBALAmaxz-GLOBALAminz)/MSDsteps;
    cout<<"Each slice will be "<<MSDint<<" Angstrom long"<<endl;
    cout <<GLOBALAminz<< " "<<GLOBALAmaxz<< " "<<(GLOBALAmaxz-GLOBALAminz)<<" "<<MSDsteps<<endl;
    MSDtics=(float*) malloc (sizeof(float)*MSDsteps*2);
    for(int i=0;i<MSDsteps;i++)
    {
      MSDtics[i*2] = GLOBALAminz +(i*MSDint);
      MSDtics[i*2+1] = GLOBALAminz +((i+1)*MSDint);
    }
    
    for(int l=0; l<(nstruct/skip);l++)
    {
    for(int i=0; i< nElt1;i++)
    {
      for(int ii=0;ii<MSDsteps;ii++)
      {
      if(A[i*3+dircover+l*nElt1*3] > MSDtics[ii*2] && A[i*3+dircover+l*nElt1*3] <= MSDtics[ii*2+1])
      {
	MSDchaos[i+l*nElt1]=ii;
      }
      }
      //cout<<i<<" "<<MSDchaos[i]<<endl;
    }
    }
   

   
  }
    
 }
 
   else if(whichmat == "box")
  {
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

    if(set_max_z == "yes")
    {
     minz = set_minz;
     maxz = set_maxz; 
    }
    if(set_max_x == "yes")
    {
     minx = set_minx;
     maxx = set_maxx;
    }
    if(set_max_y == "yes")
    {
     miny = set_miny;
     maxy = set_maxy;
    }
    cudaMalloc((void **)&dev_A,sizeof(float)*nElt1*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_exch,sizeof(int)*nElt1*(nstruct/skip));

//Copy data from host to device

    cudaMemcpy(dev_A,A,sizeof(float)*nElt1*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    

//dim3 struct to define elements of the execution configuration


    dim3 dimBlocke(32,1,32);
    dim3 dimGride((nElt1+dimBlocke.x-1)/dimBlocke.x,1,((nstruct/skip)+dimBlocke.z-1)/dimBlocke.z);

//Cuda kernal execution for distance matrix with CUDA timing API commands

    //For surface metal atoms
    cudaMemcpy(dev_exch,exch,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice);
    
    boxmat<<<dimGride,dimBlocke>>>(dev_A,dev_exch,nElt1,(nstruct/skip),maxx,minx,maxy,miny,maxz,minz);
   
    cudaMemcpy(exch,dev_exch,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyDeviceToHost);
    //For substrate atoms
     cudaFree(dev_A);
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
  } 

}
