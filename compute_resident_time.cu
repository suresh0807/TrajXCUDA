//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################


//resident time

#include "cudatools.cuh"

void compute_resident_time(void)
{

  
    
  //Aintnum is used for all cases - note
  
  //lines for exchange matrix-------------------------------------
  //------------------------------------------------------------------------------------------------------------------
cout << "metal "<< metal_species<<" counts "<<count_metal<<endl;

  int *exch1, *exch2;
  int *dev_exch1;
  exch1= (int *) malloc (sizeof(int)*nElt1*(nstruct/skip));
  exch2= (int *) malloc (sizeof(int)*nElt1*(nstruct/skip));
  Aintnum = (int *) malloc (sizeof(int)*(nstruct/skip));
  Abulknum = (int *) malloc (sizeof(int)*(nstruct/skip));
  
  
 for(int l=0; l<nstruct/skip; l++)
    {
        for(int i=0; i< nElt1; i++)
        {
          exch1[i+l*nElt1]=0;
	  exch2[i+l*nElt1]=0;
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
  float AVGintnum=0.0, AVGbulknum=0.0;
  
  int whichwater;
      if(lifetime_for=="int"){whichwater=1;}
      else if(lifetime_for=="bulk"){whichwater=1;}
      else if(lifetime_for=="all"){whichwater=2;} 
  
  
  if(lifetime_for=="int")
  {
//Allocate memory in GPU device

float *dev_A1;

    cudaMalloc((void **)&dev_B,sizeof(float)*count_metal*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_A1,sizeof(float)*nElt1*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_exch1,sizeof(int)*nElt1*(nstruct/skip));
    cudaMalloc((void **)&dev_lattice,sizeof(float)*6*(nstruct/skip));


//Copy data from host to device

    cudaMemcpy(dev_B,METAL,sizeof(float)*count_metal*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_A1,A,sizeof(float)*nElt1*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    


//dim3 struct to define elements of the execution configuration

    dim3 dimBlockc(10,10,10);
    dim3 dimGridc((nElt1+dimBlockc.x-1)/dimBlockc.x,(count_metal+dimBlockc.y-1)/dimBlockc.y,((nstruct/skip)+dimBlockc.z-1)/dimBlockc.z);

//Cuda kernal execution for distance matrix with CUDA timing API commands

  //bondist_int_bulk refers to the basic allowance within which it the correlation is allowed to start..
    cudaMemcpy(dev_exch1,exch1,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice);
    cudaMemcpy(dev_lattice,lattice,sizeof(int)*6*(nstruct/skip),cudaMemcpyHostToDevice);
    if (cell_type == "orthorhombic")
    {
    //covermat<<<dimGridc,dimBlockc>>>(dev_A1,dev_B,dev_exch1,bondist_int_bulk,nElt1,count_metal,(nstruct/skip),dev_lattice,dircover);
    covermat<<<dimGridc,dimBlockc>>>(dev_A1,dev_B,dev_exch1,bondist_int_bulk,nElt1,count_metal,(nstruct/skip),dev_lattice);
    }
    else if (cell_type == "monoclinic")
    {
     //covermatmono<<<dimGridc,dimBlockc>>>(dev_A1,dev_B,dev_exch1,bondist_int_bulk,nElt1,count_metal,(nstruct/skip),dev_lattice,dircover);
     covermatmono<<<dimGridc,dimBlockc>>>(dev_A1,dev_B,dev_exch1,bondist_int_bulk,nElt1,count_metal,(nstruct/skip),dev_lattice);
    }
    cudaMemcpy(exch1,dev_exch1,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyDeviceToHost);
  //bondist_int_bulk2 refers to the larger allowance beyong which it is definitely zero...
    cudaMemcpy(dev_exch1,exch2,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice);
    if (cell_type == "orthorhombic")
    {
    //covermat<<<dimGridc,dimBlockc>>>(dev_A1,dev_B,dev_exch1,bondist_int_bulk2,nElt1,count_metal,(nstruct/skip),dev_lattice,dircover);
    covermat<<<dimGridc,dimBlockc>>>(dev_A1,dev_B,dev_exch1,bondist_int_bulk2,nElt1,count_metal,(nstruct/skip),dev_lattice);
      
    }
    else if (cell_type == "monoclinic")
    {
     //covermatmono<<<dimGridc,dimBlockc>>>(dev_A1,dev_B,dev_exch1,bondist_int_bulk2,nElt1,count_metal,(nstruct/skip),dev_lattice,dircover);
     covermat<<<dimGridc,dimBlockc>>>(dev_A1,dev_B,dev_exch1,bondist_int_bulk2,nElt1,count_metal,(nstruct/skip),dev_lattice);
    }
    cudaMemcpy(exch2,dev_exch1,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyDeviceToHost);
    
     cudaFree(dev_A1);
     cudaFree(dev_B);
     cudaFree(dev_exch1);
  
// 1 in exch matrix means the atom satisfies the criterium
   
   for(int l=0; l<nstruct/skip; l++)
    {
      for(int i=0; i< nElt1; i++)
        {
	  if (exch2[i+l*nElt1]>0) {exch2[i+l*nElt1]=1;}
          if (exch1[i+l*nElt1]>0) {exch1[i+l*nElt1]=1; Aintnum[l]++;}
          else if (exch1[i+l*nElt1]==0) Abulknum[l]++;
        }
        AVGintnum+=Aintnum[l];
        AVGbulknum+=Abulknum[l];
    }
AVGintnum/=float(nstruct/skip);
AVGbulknum/=float(nstruct/skip);
cout <<AVGintnum<<" "<<AVGbulknum<<" "<<AVGintnum+AVGbulknum<<endl;
  }
  
  else if(lifetime_for=="bulk")
  {
 float *dev_A1;

    cudaMalloc((void **)&dev_A1,sizeof(float)*nElt1*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_exch1,sizeof(int)*nElt1*(nstruct/skip));
    cudaMalloc((void **)&dev_lattice,sizeof(float)*6*(nstruct/skip));


//Copy data from host to device

    cudaMemcpy(dev_A1,A,sizeof(float)*nElt1*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    


//dim3 struct to define elements of the execution configuration

    dim3 dimBlockc(32,1,32);
    dim3 dimGridc((nElt1+dimBlockc.x-1)/dimBlockc.x,1,((nstruct/skip)+dimBlockc.z-1)/dimBlockc.z);

//Cuda kernal execution for distance matrix with CUDA timing API commands

  
    cudaMemcpy(dev_exch1,exch1,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice);
    cudaMemcpy(dev_lattice,lattice,sizeof(int)*6*(nstruct/skip),cudaMemcpyHostToDevice);
    if (cell_type == "orthorhombic")
    {
    covermat<<<dimGridc,dimBlockc>>>(dev_A1,dev_exch1,nElt1,(nstruct/skip),MINAz,MAXAz,dircover);
    }
    else if (cell_type == "monoclinic")
    {
     covermatmono<<<dimGridc,dimBlockc>>>(dev_A1,dev_exch1,nElt1,(nstruct/skip),MINAz,MAXAz,dircover);
    }
    cudaMemcpy(exch1,dev_exch1,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyDeviceToHost);
  
    cudaMemcpy(dev_exch1,exch2,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice);
    if (cell_type == "orthorhombic")
    {
    covermat<<<dimGridc,dimBlockc>>>(dev_A1,dev_exch1,nElt1,(nstruct/skip),MIDAz-(bondist_int_bulk2/2.0),MIDAz+(bondist_int_bulk2/2.0),dircover);
    }
    else if (cell_type == "monoclinic")
    {
     covermatmono<<<dimGridc,dimBlockc>>>(dev_A1,dev_exch1,nElt1,(nstruct/skip),MIDAz-(bondist_int_bulk2/2.0),MIDAz+(bondist_int_bulk2/2.0),dircover);
    }
    cudaMemcpy(exch2,dev_exch1,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyDeviceToHost);
    
     cudaFree(dev_A1);
     cudaFree(dev_exch1);
    
// 1 in exch matrix means the atom satisfies the criterium
   
   
    for(int l=0; l<nstruct/skip; l++)
    {
      for(int i=0; i< nElt1; i++)
        {
	  if (exch2[i+l*nElt1]>0) {exch2[i+l*nElt1]=1;}
          if (exch1[i+l*nElt1]>0) {exch1[i+l*nElt1]=1; Aintnum[l]++;}
          else if (exch1[i+l*nElt1]==0) Abulknum[l]++;
        }
        AVGintnum+=Aintnum[l];
        AVGbulknum+=Abulknum[l];
    }
AVGintnum/=float(nstruct/skip);
AVGbulknum/=float(nstruct/skip);
cout <<AVGintnum<<" "<<AVGbulknum<<" "<<AVGintnum+AVGbulknum<<endl;   
  }
  
  else if (lifetime_for=="all")
  {
    cudaMalloc((void **)&dev_lattice,sizeof(float)*6*(nstruct/skip));
    cudaMemcpy(dev_lattice,lattice,sizeof(int)*6*(nstruct/skip),cudaMemcpyHostToDevice);
  }
  
  ///////////////////////////////////////////Special atom selection between 2 distances/////////////////////////////////////////////////////////////////////////

  
  
else if(choose_atoms=="yes")
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
  
 cout<<"Average surface atoms position in "<<choose_dirn<<" is "<<avgsurf<<endl;
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
    cudaMemcpy(dev_lattice,lattice,sizeof(int)*6*(nstruct/skip),cudaMemcpyHostToDevice);
    if (cell_type == "orthorhombic") 
    {
    covermat<<<dimGride,dimBlocke>>>(dev_A,dev_B,dev_exch1,minbondist,maxbondist,nElt1,count_metal,(nstruct/skip),dev_lattice,dirn,avgsurf,mintop,maxtop);
    }
    else if (cell_type == "monoclinic")
    {
    covermatmono<<<dimGride,dimBlocke>>>(dev_A,dev_B,dev_exch1,minbondist,maxbondist,nElt1,count_metal,(nstruct/skip),dev_lattice,dirn,avgsurf,mintop,maxtop);
    }
    
    
    cudaMemcpy(exch1,dev_exch1,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyDeviceToHost);
    
    cudaMemcpy(dev_exch1,exch2,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice);
    if (cell_type == "orthorhombic") 
    {
    covermat<<<dimGride,dimBlocke>>>(dev_A,dev_B,dev_exch1,minbondist,bondist_int_bulk2,nElt1,count_metal,(nstruct/skip),dev_lattice,dirn,avgsurf,mintop,maxtop);
    }
    else if (cell_type == "monoclinic")
    {
    covermatmono<<<dimGride,dimBlocke>>>(dev_A,dev_B,dev_exch1,minbondist,bondist_int_bulk2,nElt1,count_metal,(nstruct/skip),dev_lattice,dirn,avgsurf,mintop,maxtop);
    }
    
    
    cudaMemcpy(exch2,dev_exch1,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyDeviceToHost);
    
    cudaFree(dev_A);
    cudaFree(dev_B);
    cudaFree(dev_exch1);
    
    //ofstream exchout;
    //exchout.open("exchout.data");
    
    
    for(int l=0; l<nstruct/skip; l++)
      {
	//exchout<<l<<endl;
      for(int i=0; i< nElt1; i++)
        {
	  //exchout<<i<<" "<<exch[i+l*nElt1]<<endl;
	  if (exch2[i+l*nElt1]>0) {exch2[i+l*nElt1]=1;}
          if (exch1[i+l*nElt1]>0) {exch1[i+l*nElt1]=1; Aintnum[l]++;}
          else if (exch1[i+l*nElt1]==0) Abulknum[l]++;
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
  ///////////////////////////////////////////Special atom selection between 2 distances from a particlar index/////////////////////////////////////////////////////////////////////////

else if(choose_atoms_from=="yes")
{
cout << Elt1 <<" within "<< minbondist << " and " << maxbondist<<" angstrom of "<< check_index<<endl;
  
  float *index_xyz;
  index_xyz = (float *) malloc (sizeof(float)*(nstruct/skip)*3);
  
  for(int i=0 ; i<(nstruct/skip); i++)
  {
  for(int j=0 ; j<3; j++)
  {
    index_xyz[i*3+j]+=METAL[check_index*3+i*count_metal*3+j];
  }
  }
 
  
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
    cudaMalloc((void **)&dev_B,sizeof(float)*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_exch1,sizeof(int)*nElt1*(nstruct/skip));
    cudaMalloc((void **)&dev_lattice,sizeof(float)*6*(nstruct/skip));
//Copy data from host to device

    cudaMemcpy(dev_A,A,sizeof(float)*nElt1*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_B,index_xyz,sizeof(float)*(nstruct/skip)*3,cudaMemcpyHostToDevice);

//dim3 struct to define elements of the execution configuration


    dim3 dimBlocke(32,1,32);
    dim3 dimGride((nElt1+dimBlocke.x-1)/dimBlocke.x,1,((nstruct/skip)+dimBlocke.z-1)/dimBlocke.z);

//Cuda kernal execution for distance matrix with CUDA timing API commands

    //For surface metal atoms
    cudaMemcpy(dev_exch1,exch1,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice);
    cudaMemcpy(dev_lattice,lattice,sizeof(int)*6*(nstruct/skip),cudaMemcpyHostToDevice);
    if (cell_type == "orthorhombic") 
    {
      covermat<<<dimGride,dimBlocke>>>(dev_A,dev_B,dev_exch1,minbondist,maxbondist,check_index,dev_lattice,nElt1,(nstruct/skip));
    }
    else if (cell_type == "monoclinic")
    {
      covermatmono<<<dimGride,dimBlocke>>>(dev_A,dev_B,dev_exch1,minbondist,maxbondist,check_index,dev_lattice,nElt1,(nstruct/skip));
    }
    
    
    cudaMemcpy(exch1,dev_exch1,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyDeviceToHost);
    //we are using just dev_exch1 - careful- dont change the below code
    cudaMemcpy(dev_exch1,exch2,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice);
    if (cell_type == "orthorhombic") 
    {
      covermat<<<dimGride,dimBlocke>>>(dev_A,dev_B,dev_exch1,minbondist,bondist_int_bulk2,check_index,dev_lattice,nElt1,(nstruct/skip));
    }
    else if (cell_type == "monoclinic")
    {
      covermatmono<<<dimGride,dimBlocke>>>(dev_A,dev_B,dev_exch1,minbondist,bondist_int_bulk2,check_index,dev_lattice,nElt1,(nstruct/skip));
    }
    
    
    cudaMemcpy(exch2,dev_exch1,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyDeviceToHost);
    
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
	  if (exch2[i+l*nElt1]>0) {exch2[i+l*nElt1]=1;}
          if (exch1[i+l*nElt1]>0) {exch1[i+l*nElt1]=1; Aintnum[l]++;}
          else if (exch1[i+l*nElt1]==0) Abulknum[l]++;
	  //exchout<<i<<" "<<exch[i+l*nElt1]<<endl;
        }
        //exchout<<endl;
      }
      
           //exchout.close();

   
}          
   
///////////////////////////////////////////////special selection done - stored in exch///////////////////////////////////////////////////////////////////////////////////

  
////////////////////////////////////LIFETIME-PART///////////////////////////////////////////
  
  
  
  
    total_time = (nstruct/skip)*timestep;
    cout << "Each frame resolves to "<<timestep<< " fs evolution"<<endl;
    cout << "Total simulation time : " << total_time/1000000 <<" ns"<< endl;
    num_bins =  int ((nstruct/skip)/((diffuse_time*1000)/timestep));
    num_bins -= 1; 
    cout << "The analysis will be done for every "<< diffuse_time <<" ps"<<endl;
    cout <<num_bins <<" number of bins are chosen for the Hydrogen bond analysis"<<endl;
    SD_store = int((diffuse_time*1000)/timestep);
    origins = SD_store;
    cout <<"This will correlate for "<<SD_store<<" frames in the input trajectory "<<endl;
    int restarts = origins/skips;
    cout <<restarts<<" restarts from each bin will be taken"<<endl;
  
    
       SDsum1=(float *) malloc (sizeof(float)*SD_store);
        for(int ia=0;ia<SD_store;ia++)
        {
	   SDsum1[ia]=0.0;
	}
	
	float *FDsum1;
	FDsum1=(float *) malloc (sizeof(float)*SD_store);
        for(int ia=1;ia<SD_store-1;ia++)
        {
	   FDsum1[ia]=0.0;
	}
	float *Hexch1, *Hexch2;

	Hexch1 = (float *) malloc (sizeof(float)*nElt1*SD_store*2);
	Hexch2 = (float *) malloc (sizeof(float)*nElt1*SD_store*2);
 //########################################################################################################################################
//########################################################################################################################################
//BINNING AND AVERAGING THE HBAF
//#########################################################################################################################################
//######################################################################################################################################### 
float *Hexch_sized1;// correlation
float *Hexch_sized;// correlation

    for(int i=0; i<num_bins;i++)  //start bin
    {//int avgcounter=0; 
        printf("Cuda start: All lifetimes  %d \n",i);


    
         for(int i1=0;i1<SD_store*2;i1++)
         {
	     for(int k1=0; k1<nElt1; k1++)//shared hydrogen
             {
	       Hexch1[k1+i1*nElt1]=0.0;
	       Hexch2[k1+i1*nElt1]=0.0;
	     }
	 }
	 

int **initHBcol1;// storing the location of the atom in the exch matrix 

initHBcol1 = (int **) malloc (sizeof(int *)*SD_store*2);
     for(int i1=0; i1<SD_store*2; i1++)// for each frame
     {
     initHBcol1[i1] = (int *) malloc (sizeof(int)*Aintnum[i1+i*SD_store]);
     for(int j1=0; j1<Aintnum[i1+i*SD_store]; j1++)// for number of H bonds this frame
     {
      initHBcol1[i1][j1]=0;
     }
     }
     
    if(whichwater != 2)
    {
    for(int i1=0;i1<SD_store*2;i1++)
     {int chker=0;
      for(int j1=0; j1<nElt1; j1++)
       {
	if(exch1[j1+i1*nElt1+i*SD_store*nElt1] == 1 ) 
        {
	  Hexch1[j1+i1*nElt1]=1.0; //acceptor check
	  initHBcol1[i1][chker] = j1+i1*nElt1; chker+=1;
	}
	if(exch2[j1+i1*nElt1+i*SD_store*nElt1] == 1 ) 
        {
	  Hexch2[j1+i1*nElt1]=1.0; 
	}
      }
      //if(i==0){cout<<i1<<" "<<chker<<endl;}
      
     }
    }
	//cout <<"exch matrix done "<<endl;
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
     //////////////////////////////////////////////////////origin//////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
     for(int j=0; j<SD_store; j=j+skips) //go through restart points or origins
     {      
       //cout <<Aintnum[j+i*SD_store]<<" ";
       //cout<<initHBnum[j+i*SD_store]<<endl;
       if(Aintnum[j+i*SD_store] > 0){
       Hexch_sized1 = (float *) malloc (sizeof(float)*Aintnum[j+i*SD_store]*SD_store);
       Hexch_sized = (float *) malloc (sizeof(float)*Aintnum[j+i*SD_store]*SD_store);
         for(int ia=0;ia<SD_store;ia++)
         {
           for(int ja=0; ja<Aintnum[j+i*SD_store]; ja++)////only those present in the first frame in the restart bin
           { 
	        Hexch_sized1[ja+ia*Aintnum[j+i*SD_store]] = Hexch1[initHBcol1[j][ja]+ia*nElt1];
	        Hexch_sized[ja+ia*Aintnum[j+i*SD_store]] = Hexch2[initHBcol1[j][ja]+ia*nElt1];
	   }
	 }
	/* 
	 if(j==2 && i==0)
	 {
	   cout<<endl;
	   for(int ia=0;ia<SD_store;ia++)
           {
           for(int ja=0; ja<Aintnum[j+i*SD_store]; ja++)////only those present in the first frame in the restart bin
           { 
	       cout<<Hexch_sized[ja+ia*Aintnum[j+i*SD_store]] <<" ";
	   }
	   cout<<endl;
	   }
	   cout<<"after: "<<endl;
	 }
	 */
	 
	 if(HB_lifestyle=="continuous")
	 {
         for(int ja=0; ja<Aintnum[j+i*SD_store]; ja++)///for continuous lifetime
         {
           for(int ia=1;ia<SD_store;ia++)
           {
	     if(Hexch_sized[ja+(ia-1)*Aintnum[j+i*SD_store]] == 0) {Hexch_sized[ja+ia*Aintnum[j+i*SD_store]] = 0;Hexch_sized1[ja+ia*Aintnum[j+i*SD_store]] = 0;} 
	   }
	 }
	 }
	 
	 
	 if(HB_lifestyle=="transient")
	 {
	   //do transient time approximation here
	   
         for(int ja=0; ja<Aintnum[j+i*SD_store]; ja++)///for transient intermittent lifetime
         {
           for(int ia=0;ia<SD_store;ia++)
           {
	     if(Hexch_sized[ja+ia*Aintnum[j+i*SD_store]] == 0) 
	     {
	       for(int chk=1; chk<=transtime; chk++)
	       {
		 if(ia + chk >= SD_store){break;}
		 else if(Hexch_sized[ja+(ia+chk)*Aintnum[j+i*SD_store]] == 1){ Hexch_sized[ja+ia*Aintnum[j+i*SD_store]]=1;Hexch_sized1[ja+ia*Aintnum[j+i*SD_store]]=1;break;}
	       }
	    } 
	   }
	 }
	 
	 for(int ja=0; ja<Aintnum[j+i*SD_store]; ja++)///for imposing continuous lifetime after transient time
         {
           for(int ia=1;ia<SD_store;ia++)
           {
	     if(Hexch_sized[ja+(ia-1)*Aintnum[j+i*SD_store]] == 0) {Hexch_sized[ja+ia*Aintnum[j+i*SD_store]] = 0;Hexch_sized1[ja+ia*Aintnum[j+i*SD_store]] = 0;} 
	   }
	 }
	 }
      /*
       * 
       * 
       * 
         if(j==2 && i==0)
	 {
	   cout<<endl;
	   for(int ia=0;ia<SD_store;ia++)
           {
           for(int ja=0; ja<Aintnum[j+i*SD_store]; ja++)////only those present in the first frame in the restart bin
           { 
	       cout<<Hexch_sized[ja+ia*Aintnum[j+i*SD_store]] <<" ";
	   }
	   cout<<endl;
	   }
	 }
      
      
         */
       //SD=(float *) malloc (sizeof(float)*initHBnum[j+i*SD_store]*SD_store);
      
      
      
      SD1=(float *) malloc (sizeof(float)*Aintnum[j+i*SD_store]*SD_store);
        for(int ia=0;ia<SD_store;ia++)
        {
	 for(int ja=0; ja<Aintnum[j+i*SD_store]; ja++)
	 {
           //SD[ja+ia*initHBnum[j+i*SD_store]]=0.0;
	   SD1[ja+ia*Aintnum[j+i*SD_store]]=0.0;
	 }
	}
       
      
      dim3 dimBlocka(32,1,32);
      dim3 dimGrida((Aintnum[j+i*SD_store]+dimBlocka.x-1)/dimBlocka.x,1,(SD_store+dimBlocka.z-1)/dimBlocka.z);

       
      
      
      float *dev_A1;
            
      cudaMalloc((void **)&dev_A,sizeof(float)*Aintnum[j+i*SD_store]*SD_store);
      cudaMalloc((void **)&dev_A1,sizeof(float)*Aintnum[j+i*SD_store]*SD_store);
      cudaMalloc((void **)&dev_SD1,sizeof(float)*Aintnum[j+i*SD_store]*SD_store);
      
      cudaMemcpy(dev_A,Hexch_sized1,sizeof(float)*Aintnum[j+i*SD_store]*SD_store,cudaMemcpyHostToDevice);
      cudaMemcpy(dev_A1,Hexch_sized,sizeof(float)*Aintnum[j+i*SD_store]*SD_store,cudaMemcpyHostToDevice);
      cudaMemcpy(dev_SD1,SD1,sizeof(float)*Aintnum[j+i*SD_store]*SD_store,cudaMemcpyHostToDevice);

      HBAF_calc<<<dimGrida,dimBlocka>>>(dev_A,dev_A1,dev_SD1,SD_store,Aintnum[j+i*SD_store],i,j,origins,skips);

      //cudaMemcpy(SD1,dev_SD1,sizeof(float)*Aintnum[j+i*SD_store]*SD_store,cudaMemcpyDeviceToHost);
      cudaFree(dev_A);
      cudaFree(dev_A1);
      cudaFree(dev_lattice);
      //SDavg=(float *) malloc (sizeof(float)*SD_store);
      SDavg1=(float *) malloc (sizeof(float)*SD_store);
      for(int ja =0; ja<SD_store;ja++)
      {
      //SDavg[ja]=0.0;
      SDavg1[ja]=0.0;
      }
      int fairy=0;
      //cudaMalloc((void **)&dev_SDavg,sizeof(float)*SD_store);
      //cudaMemcpy(dev_SDavg,SDavg,sizeof(float)*SD_store,cudaMemcpyHostToDevice);

      dim3 dimBlocks(1024,1,1);
      dim3 dimGrids((SD_store+dimBlocks.x-1)/dimBlocks.x,1,1);
      //SDreduce<<<dimGrids,dimBlocks>>>(dev_SD,dev_SDavg,SD_store,initHBnum[j+i*SD_store],fairy);
      //cudaMemcpy(SDavg,dev_SDavg,sizeof(float)*SD_store,cudaMemcpyDeviceToHost);
       
      //cudaFree(dev_SD);
      //cudaFree(dev_SDavg);
      
      cudaMalloc((void **)&dev_SDavg1,sizeof(float)*SD_store);
      cudaMemcpy(dev_SDavg1,SDavg1,sizeof(float)*SD_store,cudaMemcpyHostToDevice);
      SDreduce<<<dimGrids,dimBlocks>>>(dev_SD1,dev_SDavg1,SD_store,Aintnum[j+i*SD_store],fairy);
      cudaMemcpy(SDavg1,dev_SDavg1,sizeof(float)*SD_store,cudaMemcpyDeviceToHost);
       
      cudaFree(dev_SD1);
      cudaFree(dev_SDavg1);
      
           
     // SDsum1[0]+=1;
//cout<<Aintnum[j+i*SD_store]<<" "<<SDavg1[0]<<" ";
      for(int ja =0; ja<SD_store;ja++)
      {
      SDsum1[ja]+=SDavg1[ja];
      //SDsum1[ja]+=(SDavg1[ja]/SDavg1[0]);
      }
  
      for(int ja=1; ja<SD_store -1;ja++)
      {
      //FDsum[ja]+=-((SDavg[ja+1]-SDavg[ja-1]) / (((ja+1)*(timestep/1000))-((ja-1)*(timestep/1000))));
      FDsum1[ja]+=-((SDavg1[ja+1]-SDavg1[ja-1]) / (((ja+1)*(timestep/1000))-((ja-1)*(timestep/1000))));      
      }
      
      
      
      
      free(Hexch_sized);
      //free(SDavg);
      //free(SD);
      free(Hexch_sized1);
      free(SDavg1);
      //free(SD1);
      //avgcounter++;
       }
    }//origins over
    

  
  free(initHBcol1);
     
   }

  free(Hexch1);
   free(Hexch2);
  ofstream Hexchplot, Hexchplotac;
  
      Hexchplotac.open("ct.data");
     for(int ja =0; ja<SD_store;ja++)
     {
       //SDsum[ja]/=float(num_bins*restarts);
       SDsum1[ja]/=float(num_bins*restarts);
       //Hexchplotac <<ja*timestep/1000<<" "<<SDsum1[ja]<<endl;
       Hexchplotac <<ja*timestep/1000<<" "<<SDsum1[ja]/SDsum1[0]<<endl;
       //Hexchplot <<ja*timestep/1000<<" "<<SDsum[ja]<<endl;
      } 
      //Hexchplot.close();
      Hexchplotac.close();
      
    
      //Hexchplot.open("ft-pair.data");
      Hexchplotac.open("ft.data");
     for(int ja =1; ja<SD_store-1;ja++)
      {
       //FDsum[ja]/=num_bins*restarts;
       FDsum1[ja]/=num_bins*restarts;
       //Hexchplot <<ja*timestep/1000<<" "<<FDsum[ja]<<endl;
       Hexchplotac <<ja*timestep/1000<<" "<<FDsum1[ja]<<endl;
      }  
      //Hexchplot.close();
      Hexchplotac.close();
     
     int avg_every = 100;
     float *FD_avg1;
     //FD_avg=(float *) malloc (sizeof(float)*SD_store);
     FD_avg1=(float *) malloc (sizeof(float)*SD_store);
        for(int ia=0;ia<SD_store;ia++)
        {
          // FD_avg[ia]=0.0;
	   FD_avg1[ia]=0.0;
	}
     
      for(int ja=(avg_every/2)+1; ja<SD_store -(avg_every/2)-1;ja++)
      {
	for(int ka=ja-(avg_every/2);ka<=ja+(avg_every/2)-1;ka++)
	{
	 // FD_avg[ja]+= FDsum[ka];
	  FD_avg1[ja]+= FDsum1[ka];
	}
	//FD_avg[ja]/=avg_every;
	FD_avg1[ja]/=avg_every;
      }
      
      Hexchplotac.open("ft-avg.data");
      //Hexchplot.open("ft-avg-pair.data");
     for(int ja =(avg_every/2)+1; ja<SD_store -(avg_every/2)-1;ja++)
      {
       //Hexchplot <<ja*timestep/1000<<" "<<FD_avg[ja]<<endl;
       Hexchplotac <<ja*timestep/1000<<" "<<FD_avg1[ja]<<endl;
      }  
      //Hexchplot.close();
      Hexchplotac.close();
     
      
      
      
     float lifetime1;
     float *int_SDsum1;
     float *cum_SDsum1;
     //int_SDsum = (float*) malloc (sizeof(float)*SD_store);
     //cum_SDsum = (float*) malloc (sizeof(float)*SD_store);
     int_SDsum1 = (float*) malloc (sizeof(float)*SD_store);
     cum_SDsum1 = (float*) malloc (sizeof(float)*SD_store);
        for(int ia=0;ia<SD_store;ia++)
        {
           //int_SDsum[ia]=0.0;
	   //cum_SDsum[ia]=0.0;
	   int_SDsum1[ia]=0.0;
	   cum_SDsum1[ia]=0.0;
	}
     for(int ja =1; ja<SD_store;ja++)
      {
       //int_SDsum[ja]= (((ja*timestep*0.001) - ((ja-1)*timestep*0.001)) * ((SDsum[ja] + SDsum[ja-1])/2.0));
       //cum_SDsum[ja]=cum_SDsum[ja-1]+int_SDsum[ja];
       int_SDsum1[ja]= (((ja*timestep*0.001) - ((ja-1)*timestep*0.001)) * (((SDsum1[ja]/SDsum1[0]) + (SDsum1[ja-1]/SDsum1[0]))/2.0));
       cum_SDsum1[ja]=cum_SDsum1[ja-1]+int_SDsum1[ja];
       if(ja == SD_store-1) {//lifetime = cum_SDsum[ja];
	 lifetime1 = cum_SDsum1[ja];}
      }  
     
     //Hexchplot.open("ct-integrate-pair.data");
     Hexchplotac.open("ct-integrate.data");
     for(int ja =0; ja<SD_store;ja++)
      {
       //Hexchplot <<ja*timestep/1000<<" "<<cum_SDsum[ja]<<endl;
       Hexchplotac <<ja*timestep/1000<<" "<<cum_SDsum1[ja]<<endl;
      }  
      //Hexchplot.close();
      Hexchplotac.close();
     //Hexchplot.open("Hbond-lifetime-pair.data");
     //Hexchplot <<"Lifetime from the integral of c(t) is : "<<lifetime<<" ps"<<endl;
     //Hexchplot.close();
     Hexchplot.open("lifetime.data");
     Hexchplot <<"Lifetime from the integral of c(t) is : "<<lifetime1<<" ps"<<endl;
     Hexchplot.close();
     
 

 
  free(A);
  free(B);
 }
