//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################


///Now can process a very long trajectory also!! Still some optimizations needed in the initial regions

#include "cudatools.cuh"

void compute_adatoms(void)
{

  //lines for exchange matrix------NO implication whatsoever for the Hbonds computation-------------------------------
  //------------------------------------------------------------------------------------------------------------------
cout << "metal "<< metal_species<<" counts "<<count_metal<<endl;

  int *exch1;
  int *dev_exch1;
  exch1= (int *) malloc (sizeof(int)*nElt1*(nstruct/skip));
  Aintnum = (int *) malloc (sizeof(int)*(nstruct/skip));
  Abulknum = (int *) malloc (sizeof(int)*(nstruct/skip));
  
  // Initialize distance matrix and histogram matrix
 for(int l=0; l<nstruct/skip; l++)
    {
        for(int i=0; i< nElt1; i++)
        {
          exch1[i+l*nElt1]=0;
        }
        Aintnum[l]=0;
        Abulknum[l]=0;
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

  MAXAz=MIDAz+(bondist_int_bulk/2.0);
  MINAz=MIDAz-(bondist_int_bulk/2.0);
  
 
   ///////////////////////////////////////////Special atom selection between 2 distances/////////////////////////////////////////////////////////////////////////

  
if(choose_atoms=="yes")
{
  float surf,avgsurf=0.0;
cout << Elt1 <<" within "<< minbondist << " and " << maxbondist<<" angstrom of "<< metal_species<<endl;
  

  int dirn;
  if(choose_dirn=="z") dirn=2;
  else if(choose_dirn=="x") dirn=0;
  else if(choose_dirn=="y") dirn=1;
  
  
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
  
 cout<<"Average surface atom position in "<<choose_dirn<<" is "<<avgsurf<<endl;
 cout<<"Average number of surface atoms is "<< surfatom<<endl;
  
 
//Allocate memory in GPU device

    
    cudaMalloc((void **)&dev_A,sizeof(float)*nElt1*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_B,sizeof(float)*count_metal*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_exch1,sizeof(int)*nElt1*(nstruct/skip));
    cudaMalloc((void **)&dev_lattice,sizeof(float)*6*(nstruct/skip));
//Copy data from host to device

    cudaMemcpy(dev_A,A,sizeof(float)*nElt1*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_B,METAL,sizeof(float)*count_metal*(nstruct/skip)*3,cudaMemcpyHostToDevice);

//dim3 struct to define elements of the execution configuration


    dim3 dimBlocke(10,10,10);
    dim3 dimGride((nElt1+dimBlocke.x-1)/dimBlocke.x,(count_metal+dimBlocke.y-1)/dimBlocke.y,((nstruct/skip)+dimBlocke.z-1)/dimBlocke.z);

//Cuda kernal execution for distance matrix with CUDA timing API commands

    //For surface metal atoms
    cudaMemcpy(dev_exch1,exch1,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice);

    if (cell_type == "orthorhombic") 
    {
    covermat<<<dimGride,dimBlocke>>>(dev_A,dev_B,dev_exch1,minbondist,maxbondist,nElt1,count_metal,(nstruct/skip),dev_lattice,mintop,maxtop,dirn,avgsurf);
    }
    else if (cell_type == "monoclinic")
    {
    covermatmono<<<dimGride,dimBlocke>>>(dev_A,dev_B,dev_exch1,minbondist,maxbondist,nElt1,count_metal,(nstruct/skip),dev_lattice,mintop,maxtop,dirn,avgsurf);
    }
    
    
    cudaMemcpy(exch1,dev_exch1,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyDeviceToHost);
    
    
    cudaFree(dev_A);
    cudaFree(dev_B);
    
    //ofstream exchout;
    //exchout.open("exchout.data");
    
    
    for(int l=0; l<nstruct/skip; l++)
      {
	//exchout<<l<<endl;
      for(int i=0; i< nElt1; i++)
        {
	  //exchout<<i<<" "<<exch[i+l*nElt1]<<endl;
          if (exch1[i+l*nElt1]>0) {exch1[i+l*nElt1]=1;}
	  //exchout<<i<<" "<<exch[i+l*nElt1]<<endl;
        }
        //exchout<<endl;
      }
      
      //exchout.close();
    // cudaMemcpy(dev_exch1,exch1,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice); 
     
  
}          
   
///////////////////////////////////////////////special selection done - stored in exch///////////////////////////////////////////////////////////////////////////////////

//exchange matrix created--------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------
 /* debug--------------
 for(int l=0; l<nstruct/skip; l++)
    {
        for(int i=0; i< nElt1; i++)
        {
          if(exch1[i+l*nElt1]==1){cout<<l<<" "<<i<<endl;}
        }
    }
  */
  int *Onum; //number of O atoms within a specific distance (must be greater than 2 if used 3.2) of the H atom.
  int *Hbondneigh; // Ids of the (O) neighbours of the H atom. 
  float *Hbondneighdist; // OH distance for all neighbours.
  
  float *Hdistmat, *dev_Hdistmat;// Distance matrix to capture the O within x \AA of H atom.
  
  int OHsamples=50; /// maximum number of OH bonds considered within bonddist sphere
  
    Onum = (int *) malloc (sizeof(int)*nElt1*(nstruct/skip));
    Hbondneigh = (int *) malloc (sizeof(int)*OHsamples*nElt1*(nstruct/skip));
    Hbondneighdist = (float *) malloc (sizeof(float)*OHsamples*nElt1*(nstruct/skip));
    Hdistmat = (float *) malloc (sizeof(float)*nElt1*nElt2);
    
  for(int i=0;i<nstruct/skip;i++)
  {
    for(int j=0; j<nElt1; j++)
    {
      Onum[j+i*nElt1] =0;
      for(int k =0 ;k<OHsamples;k++)
      {
      Hbondneigh[j*OHsamples+i*nElt1*OHsamples+k]=0;
      Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k]=0.0;
      }
    }
  }   
      
 for(int i=0;i<nElt1;i++)
  {
    for(int j=0; j<nElt2; j++)
    {
      Hdistmat[j+i*nElt2]=0.0;
    }
  }
      
cout <<"gpu begins distance computation"<<endl;
//Allocate memory in GPU device

    cudaMalloc((void **)&dev_A,sizeof(float)*nElt1*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_B,sizeof(float)*nElt2*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_Hdistmat,sizeof(float)*nElt1*nElt2);
    
//Copy data from host to device

    cudaMemcpy(dev_A,A,sizeof(float)*nElt1*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_B,B,sizeof(float)*nElt2*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_Hdistmat,Hdistmat,sizeof(float)*nElt1*nElt2,cudaMemcpyHostToDevice);
    
   //dim3 struct to define elements of the execution configuration


    dim3 dimBlock(32,32,1);
    dim3 dimGrid((nElt1+dimBlock.x-1)/dimBlock.x,(nElt2+dimBlock.y-1)/dimBlock.y,1);
for(int i = 0; i<(nstruct/skip); i++) /// go into each frame
{
   for(int i1=0;i1<nElt1;i1++)///initialize Hdistmat for 1 frame
  {
    for(int j1=0; j1<nElt2; j1++)
    {
      Hdistmat[j1+i1*nElt2]=0.0;
    }
  }
  cudaMemcpy(dev_Hdistmat,Hdistmat,sizeof(float)*nElt1*nElt2,cudaMemcpyHostToDevice);
//Cuda kernal execution for distance matrix with CUDA timing API commands
    if(cell_type == "orthorhombic")
    {
          Hbondmat<<<dimGrid,dimBlock>>>(dev_A,dev_B,dev_exch1,dev_Hdistmat,bondist,nElt1,nElt2,dev_lattice,i);
    }
    else if(cell_type == "monoclinic")
    {
          Hbondmatmono<<<dimGrid,dimBlock>>>(dev_A,dev_B,dev_exch1,dev_Hdistmat,bondist,nElt1,nElt2,dev_lattice,i);      
    }

cudaMemcpy(Hdistmat,dev_Hdistmat,sizeof(float)*nElt1*nElt2,cudaMemcpyDeviceToHost);
   
  for(int i1=0;i1<nElt1;i1++) //hydrogen
  {
    int k=0;
    for(int j1=0; j1<nElt2; j1++) //oxygen 
    {
      if(Hdistmat[j1+i1*nElt2] !=0.0) // check if an oxygen atom is within the 'bonddist' distance from hydrogen atom
      {
	Onum[i1+i*nElt1]++;
	Hbondneigh[i1*OHsamples+i*nElt1*OHsamples+k] = j1;
	Hbondneighdist[i1*OHsamples+i*nElt1*OHsamples+k] = Hdistmat[j1+i1*nElt2];
	k++;
      }
    }
  }

}   
    
    cudaFree(dev_A);
    cudaFree(dev_B);
    cudaFree(dev_Hdistmat);
    cudaFree(dev_lattice);
    cudaFree(dev_exch1);


//Hbondneigh contain the IDs of four possible oxygen neighbour atoms if the hydrogen chosen is participating in a hydrogen bond.
//Onum has the number of oxygen atoms within a specified distance (see input file).

cout <<"gpu ends"<<endl; 

      
//swapping the indices of nearest oxygens to get the two shortest OH bonds up front, 
//since it usually describes the H bond

cout<<"Swapping of OH bonds based on distance "<<endl;

float swapa,swapb;
    for(int i=0;i<(nstruct/skip);i++)
     {
      for(int j=0; j<nElt1; j++)
       {
	 
	for(int swf=0; swf<Onum[j+i*nElt1]-1;swf++)
	{
	for(int sws=swf+1; sws<Onum[j+i*nElt1];sws++)
	{
	if(Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+swf] > Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+sws] && sws !=swf && Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+sws] !=0.0 && Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+swf] !=0.0) 
	{
	   swapa = Hbondneigh[j*OHsamples+i*nElt1*OHsamples+swf]; 
	    Hbondneigh[j*OHsamples+i*nElt1*OHsamples+swf]=Hbondneigh[j*OHsamples+i*nElt1*OHsamples+sws];
	    Hbondneigh[j*OHsamples+i*nElt1*OHsamples+sws]=swapa;
	   swapb = Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+swf]; 
	    Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+swf]=Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+sws];
	    Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+sws]=swapb;
	}
	}
	}	
       }
     }	
   
cout<<" done !!!"<<endl;
   
   
   ofstream Hneigh;
Hneigh.open("Adneigh.data");
Hneigh<<"S H N F FD      S   SD      OHO     ALEX<   MATTI   OO"<<endl;
  for(int i=0;i<nstruct/skip;i++)
  {      
    for(int j=0; j<nElt1; j++)
    {  
      for(int k=0; k<OHsamples;k++)
      {
	if(Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k] !=0 && Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k] != Hbondneighdist[j*OHsamples+i*nElt1*OHsamples])
	{
      Hneigh<<i<<" "<<j<<" "<<Onum[j+i*nElt1]<<" "<<Hbondneigh[j*OHsamples+i*nElt1*OHsamples]<<" "<<Hbondneighdist[j*OHsamples+i*nElt1*OHsamples]<<" "<<Hbondneigh[j*OHsamples+i*nElt1*OHsamples+k]<<" "<<Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k]<<endl;
	}
      //Hneigh<<j<<" "<<Onum[j+i*nElt1]<<" "<<Hbondneigh[j*OHsamples+i*nElt1*OHsamples]<<" "<<Hbondneighdist[j*OHsamples+i*nElt1*OHsamples]<<" "<<Hbondneigh[j*OHsamples+i*nElt1*OHsamples+1]<<" "<<Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+1]<<" "<<Hbondneigh[j*OHsamples+i*nElt1*OHsamples+2]<<" "<<Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+2]<<" "<<Hbondneigh[j*OHsamples+i*nElt1*OHsamples+3]<<" "<<Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+3]<<endl;
      }
    } 
  }  
Hneigh.close();
  
cout<<"Done!!!"<<endl;

int *Batop, *Bbridge, *Bhcp, *Bfcc;
int *Catop, *Cbridge, *Chcp, *Cfcc;

Batop= (int *) malloc (sizeof(int)*nstruct/skip*nElt1);
Bbridge= (int *) malloc (sizeof(int)*nstruct/skip*nElt1);
Bhcp= (int *) malloc (sizeof(int)*nstruct/skip*nElt1);
Bfcc= (int *) malloc (sizeof(int)*nstruct/skip*nElt1);
Catop= (int *) malloc (sizeof(int)*nstruct/skip);
Cbridge= (int *) malloc (sizeof(int)*nstruct/skip);
Chcp= (int *) malloc (sizeof(int)*nstruct/skip);
Cfcc= (int *) malloc (sizeof(int)*nstruct/skip);
for(int i=0;i<nstruct/skip;i++)
  {      
    for(int j=0; j<nElt1; j++)
    {
      Batop[j+i*nElt1]=0;
      Bbridge[j+i*nElt1]=0;
      Bhcp[j+i*nElt1]=0;
      Bfcc[j+i*nElt1]=0;
    }
      Catop[i]=0;
      Cbridge[i]=0;
      Chcp[i]=0;
      Cfcc[i]=0;
  }

int atop=0;
int bridge=0;
int fcc=0;
int hcp=0;
int k=0;
for(int i=0;i<nstruct/skip;i++)
//for(int i=0;i<1;i++)
  {      
    for(int j=0; j<nElt1; j++)
    {  
      float A1=Hbondneighdist[j*OHsamples+i*nElt1*OHsamples];
      float A2=Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+1];
      float A3=Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+2];
      
      if(A1<2.75 && A2<2.75 && A3>2.70&& A1!=0&&A2!=0&&A3!=0) {bridge++; Bbridge[j+i*nElt1]=1;Cbridge[i]++;
	//cout<<"bridge "<<i<<endl;
	cout<<i<<" "<<A1<<" "<<A2<<" "<<A3<<" bridge"<<endl;
	      }
      else if(A2>2.75 && A1!=0 && A2!=0) 
        {
	atop++; Batop[j+i*nElt1]=1;Catop[i]++;
	cout<<i<<" "<<A1<<" "<<A2<<" "<<A3<<" atop"<<endl;
	}
      else if(A1<2.75 && A2<2.75 && A3<2.75&& A1!=0&&A2!=0&&A3!=0) 
      {
	
	for(k=3;k<25;k++)
	{
	  
	  if(Hbondneighdist[k+j*OHsamples+i*nElt1*OHsamples] !=0)
	  {
	    //cout<<A[j*3+i*nElt1*3]<<" "<<A[Hbondneigh[k+j*OHsamples+i*nElt1*OHsamples]*3+i*nElt1*3]<<" "<<A[j*3+i*nElt1*3+1]<<" "<<A[Hbondneigh[k+j*OHsamples+i*nElt1*OHsamples]*3+i*nElt1*3+1]<<endl;
	  if(fabs(A[j*3+i*nElt1*3]-A[Hbondneigh[k+j*OHsamples+i*nElt1*OHsamples]*3+i*nElt1*3])<0.6 && fabs(A[j*3+i*nElt1*3+1]-A[Hbondneigh[k+j*OHsamples+i*nElt1*OHsamples]*3+i*nElt1*3+1])<0.6)
	  {
	    hcp++;Bhcp[j+i*nElt1]=1;Chcp[i]++;cout<<i<<" "<<A1<<" "<<A2<<" "<<A3<<" hcp"<<endl;break;
	  }
	  }
        }
        if(k==25){fcc++; Bfcc[j+i*nElt1]=1;Cfcc[i]++;cout<<i<<" "<<A1<<" "<<A2<<" "<<A3<<" fcc"<<endl;}
      }
      else if(A1!=0&&A2!=0&&A3!=0) 
      {
	cout<<i<<endl;
      }
      //if(exch1[j+i*nElt1]>0) {cout<<i<<endl;}
    }
  }

   cout<<"Atop :"<<atop<<" Bridge :"<<bridge<<" fcc :"<<fcc<<" hcp :"<<hcp<<endl;
  free(exch1);
  free(Hbondneighdist);
  free(Onum);
  free(Hbondneigh);
  free(Hdistmat);
  
  
  compute_lifetime(Bfcc,Cfcc,25);
  system("mv ct.data ct-fcc.data; mv ct-integrate.data ct-fcc-integrate.data; mv ft.data ft-fcc.data; mv ft-avg.data ft-fcc-avg.data; mv lifetime.data lifetime-fcc.data");
  compute_lifetime(Bhcp,Chcp,25);
  system("mv ct.data ct-hcp.data; mv ct-integrate.data ct-hcp-integrate.data; mv ft.data ft-hcp.data; mv ft-avg.data ft-hcp-avg.data; mv lifetime.data lifetime-hcp.data");
  compute_lifetime(Bbridge,Cbridge,5);
  system("mv ct.data ct-bridge.data; mv ct-integrate.data ct-bridge-integrate.data; mv ft.data ft-bridge.data; mv ft-avg.data ft-bridge-avg.data; mv lifetime.data lifetime-bridge.data");
  compute_lifetime(Batop,Catop,5);
  system("mv ct.data ct-atop.data; mv ct-integrate.data ct-atop-integrate.data; mv ft.data ft-atop.data; mv ft-avg.data ft-atop-avg.data; mv lifetime.data lifetime-atop.data");
  
  free(Batop);
  free(Bbridge);
  free(Bfcc);
  free(Bhcp);
  free(Catop);
  free(Cbridge);
  free(Cfcc);
  free(Chcp);
}
