//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"

void compute_Hbonds(void)
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
  
  MAXAz=MIDAz+(bondist/2.0);
  MINAz=MIDAz-(bondist/2.0);
  float AVGintnum=0.0, AVGbulknum=0.0;
  
  if(HB_for=="int")
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

  
    cudaMemcpy(dev_exch1,exch1,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice);
    cudaMemcpy(dev_lattice,lattice,sizeof(int)*6*(nstruct/skip),cudaMemcpyHostToDevice);
    if (cell_type == "orthorhombic")
    {
    covermat<<<dimGridc,dimBlockc>>>(dev_A1,dev_B,dev_exch1,bondist_int_bulk,nElt1,count_metal,(nstruct/skip),dev_lattice);
    }
    else if (cell_type == "monoclinic")
    {
     covermatmono<<<dimGridc,dimBlockc>>>(dev_A1,dev_B,dev_exch1,bondist_int_bulk,nElt1,count_metal,(nstruct/skip),dev_lattice);
    }
    cudaMemcpy(exch1,dev_exch1,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyDeviceToHost);
  
     cudaFree(dev_A1);
     cudaFree(dev_B);
     
   /*  ofstream exchplot;
    exchplot.open("exchange1.data");

   for(int l=0; l<(nstruct/skip); l++)
    {
      exchplot <<l<<" ";
      for(int j=0; j<nElt1; j++)
      {

          exchplot << exch1[j+l*nElt1]<<" ";
      }
      exchplot <<endl;
    }
    exchplot.close();*/

   
   // 1 in exch matrix meand a HBOND is present in interface, 0 means it is in bulk
   
    for(int l=0; l<nstruct/skip; l++)
    {
      for(int i=0; i< nElt1; i++)
        {
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
  
  else if(HB_for=="bulk")
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
  
     cudaFree(dev_A1);
     
   /*  ofstream exchplot;
    exchplot.open("exchange1.data");

   for(int l=0; l<(nstruct/skip); l++)
    {
      exchplot <<l<<" ";
      for(int j=0; j<nElt1; j++)
      {

          exchplot << exch1[j+l*nElt1]<<" ";
      }
      exchplot <<endl;
    }
    exchplot.close();*/

   
   // 1 in exch matrix meand a HBOND is present in interface, 0 means it is in bulk
   
   
    for(int l=0; l<nstruct/skip; l++)
    {
      for(int i=0; i< nElt1; i++)
        {
          if (exch1[i+l*nElt1]>0) {exch1[i+l*nElt1]=1; Abulknum[l]++;}
          else if (exch1[i+l*nElt1]==0) Aintnum[l]++;
        }
        AVGintnum+=Aintnum[l];
        AVGbulknum+=Abulknum[l];
    }
AVGintnum/=float(nstruct/skip);
AVGbulknum/=float(nstruct/skip);
cout <<AVGintnum<<" "<<AVGbulknum<<" "<<AVGintnum+AVGbulknum<<endl;   
  }
  
  else if (HB_for=="all")
  {
    cudaMalloc((void **)&dev_lattice,sizeof(float)*6*(nstruct/skip));
    cudaMemcpy(dev_lattice,lattice,sizeof(int)*6*(nstruct/skip),cudaMemcpyHostToDevice);
  }
//exchange matrix created--------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------
  
  
  
  int *Onum; //number of O atoms within a specific distance (must be greater than 2 if used 3.2) of the H atom.
  int *Hbondneigh; // Ids of the (O) neighbours of the H atom. 
  float *Hbondneighdist; // OH distance for all neighbours.
  float *Hbondxyz; //storage.
  
  float *Hdistmat, *dev_Hdistmat;// Distance matrix to capture the O within x \AA of H atom.
  
  int OHsamples=20; /// maximum number of OH bonds considered within bonddist sphere
  
    Onum = (int *) malloc (sizeof(int)*nElt1*(nstruct/skip));
    Hbondneigh = (int *) malloc (sizeof(int)*OHsamples*nElt1*(nstruct/skip));
    Hbondneighdist = (float *) malloc (sizeof(float)*OHsamples*nElt1*(nstruct/skip));
    Hbondxyz = (float *) malloc (sizeof(float)*9*nElt1*(nstruct/skip));  //X,Y and Z of Hatom followed by angle (depends on criterium set), OO distance, neighbour O id, and O-H distance,OH distance,O-HO/OO-H
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
      for(int k =0 ;k<9;k++)
      {
      Hbondxyz[j*9+i*nElt1*9+k]=0.0;
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
      
cout <<"gpu begins OH distance computation"<<endl;
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
          Hbondmat<<<dimGrid,dimBlock>>>(dev_A,dev_B,dev_Hdistmat,bondist,nElt1,nElt2,dev_lattice,i);
    }
    else if(cell_type == "monoclinic")
    {
          Hbondmatmono<<<dimGrid,dimBlock>>>(dev_A,dev_B,dev_Hdistmat,bondist,nElt1,nElt2,dev_lattice,i);      
    }

cudaMemcpy(Hdistmat,dev_Hdistmat,sizeof(float)*nElt1*nElt2,cudaMemcpyDeviceToHost);
   
  for(int i1=0;i1<nElt1;i1++) //hydrogen
  {
    int k=0;
    for(int j1=0; j1<nElt2; j1++) //oxygen 
    {
      if(Hdistmat[j1+i1*nElt2] !=0.0) // check if an axygen atom is within the 'bonddist' distance from hydrogen atom
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
   
   
cout<<"All angles and distance computation begins"<<endl;
   
 float  *distOOsamples, *anglesOHOsamples, *anglesHOOsamples;
 
 distOOsamples = (float*) malloc (sizeof(float)*(OHsamples)*nElt1*(nstruct/skip));
 anglesOHOsamples = (float*) malloc (sizeof(float)*(OHsamples)*nElt1*(nstruct/skip));
 anglesHOOsamples = (float*) malloc (sizeof(float)*(OHsamples)*nElt1*(nstruct/skip)); 
 
 for(int i=0; i<(nstruct/skip); i++)
 {
   for(int j=0; j<nElt1; j++)
   {
     for(int k=0; k<(OHsamples); k++)
     {
       distOOsamples[k+j*(OHsamples)+i*nElt1*(OHsamples)]=0.0;
       anglesOHOsamples[k+j*(OHsamples)+i*nElt1*(OHsamples)]=0.0;
       anglesHOOsamples[k+j*(OHsamples)+i*nElt1*(OHsamples)]=0.0;
     }
   }
 }
 float vec1[3],vec2[3];
 for(int i=0;i<(nstruct/skip);i++) // into each frame one at a time
     {
      for(int j=0; j<nElt1; j++) // for every hydrogen atom
       {
	 for(int k =0 ;k<3;k++)
         {
	 vec1[k] = B[i*nElt2*3+Hbondneigh[j*OHsamples+i*nElt1*OHsamples]*3+k];
	 }
	 
	 //for(int l=0; l<OHsamples; l++)
	 for(int l=0; l<8; l++)//only the first eight O neighbours enough
	 {
	   
	 if(cell_type == "orthorhombic")
	 {
	for(int k =0 ;k<3;k++)
	{
	vec2[k] = B[i*nElt2*3+Hbondneigh[j*OHsamples+i*nElt1*OHsamples+l]*3+k];
	if(k==0) 
	{
	  if(fabs(A[j*3+i*nElt1*3+k]-vec1[k]) > xvec/2) 
	          {if (A[j*3+i*nElt1*3+k]-vec1[k] > 0) {vec1[k] = vec1[k] +xvec;}
		   else {vec1[k] = vec1[k] -xvec;}}
	  if(fabs(A[j*3+i*nElt1*3+k]-vec2[k]) > xvec/2) 
	          {if (A[j*3+i*nElt1*3+k]-vec2[k] > 0) {vec2[k] = vec2[k] +xvec;}
		   else {vec2[k] = vec2[k] -xvec;}}	   
	}
		   
	else if(k==1) 
	{
	  if(fabs(A[j*3+i*nElt1*3+k]-vec1[k]) > yvec/2) 
	          {if (A[j*3+i*nElt1*3+k]-vec1[k] > 0) {vec1[k] = vec1[k] +yvec;}
		   else {vec1[k] = vec1[k] -yvec;}}
	  if(fabs(A[j*3+i*nElt1*3+k]-vec2[k]) > yvec/2) 
	          {if (A[j*3+i*nElt1*3+k]-vec2[k] > 0) {vec2[k] = vec2[k] +yvec;}
		   else {vec2[k] = vec2[k] -yvec;}}   
	}
		   
	else if(k==2) 
	{
	  if(fabs(A[j*3+i*nElt1*3+k]-vec1[k]) > zvec/2) 
	          {if (A[j*3+i*nElt1*3+k]-vec1[k] > 0) {vec1[k] = vec1[k] +zvec;}
		   else {vec1[k] = vec1[k] -zvec;}}
	  if(fabs(A[j*3+i*nElt1*3+k]-vec2[k]) > zvec/2) 
	          {if (A[j*3+i*nElt1*3+k]-vec2[k] > 0) {vec2[k] = vec2[k] +zvec;}
		   else {vec2[k] = vec2[k] -zvec;}}	   
	}
	}
	 }
	 else if(cell_type =="monoclinic")
	 {
	int k=0;
	vec2[0] = B[i*nElt2*3+Hbondneigh[j*OHsamples+i*nElt1*OHsamples+l]*3];
	vec2[1] = B[i*nElt2*3+Hbondneigh[j*OHsamples+i*nElt1*OHsamples+l]*3+1];
	vec2[2] = B[i*nElt2*3+Hbondneigh[j*OHsamples+i*nElt1*OHsamples+l]*3+2];
	k=1;
	if(fabs(A[j*3+i*nElt1*3+k]-vec1[k]) > yvec/2) 
	          {if (A[j*3+i*nElt1*3+k]-vec1[k] > 0) {vec1[k] = vec1[k] +yvec; vec1[0]=vec1[0]+yxvec;}
		   else {vec1[k] = vec1[k] -yvec; vec1[0]=vec1[0]-yxvec;}}
	if(fabs(A[j*3+i*nElt1*3+k]-vec2[k]) > yvec/2) 
	          {if (A[j*3+i*nElt1*3+k]-vec2[k] > 0) {vec2[k] = vec2[k] +yvec; vec2[0]=vec2[0]+yxvec;}
		   else {vec2[k] = vec2[k] -yvec; vec2[0]=vec2[0]-yxvec;}}	   
	k=0;
	if(fabs(A[j*3+i*nElt1*3+k]-vec1[k]) > xvec/2) 
	          {if (A[j*3+i*nElt1*3+k]-vec1[k] > 0) {vec1[k] = vec1[k] +xvec;}
		   else {vec1[k] = vec1[k] -xvec;}}
	if(fabs(A[j*3+i*nElt1*3+k]-vec2[k]) > xvec/2) 
	          {if (A[j*3+i*nElt1*3+k]-vec2[k] > 0) {vec2[k] = vec2[k] +xvec;}
		   else {vec2[k] = vec2[k] -xvec;}} 	 
	k=2;
	if(fabs(A[j*3+i*nElt1*3+k]-vec1[k]) > zvec/2) 
	          {if (A[j*3+i*nElt1*3+k]-vec1[k] > 0) {vec1[k] = vec1[k] +zvec;}
		   else {vec1[k] = vec1[k] -zvec;}}
	if(fabs(A[j*3+i*nElt1*3+k]-vec2[k]) > zvec/2) 
	          {if (A[j*3+i*nElt1*3+k]-vec2[k] > 0) {vec2[k] = vec2[k] +zvec;}
		   else {vec2[k] = vec2[k] -zvec;}}
	
	 }
	 //angle computed with the central element at first and the others following
	 anglesOHOsamples[l+j*(OHsamples)+i*nElt1*(OHsamples)]= angle(A[j*3+i*nElt1*3],A[j*3+i*nElt1*3+1],A[j*3+i*nElt1*3+2],vec1[0], vec1[1], vec1[2], vec2[0], vec2[1], vec2[2]);
	 anglesHOOsamples[l+j*(OHsamples)+i*nElt1*(OHsamples)]= angle(vec1[0], vec1[1], vec1[2],vec2[0], vec2[1], vec2[2],A[j*3+i*nElt1*3],A[j*3+i*nElt1*3+1],A[j*3+i*nElt1*3+2]);
	 distOOsamples[l+j*(OHsamples)+i*nElt1*(OHsamples)]= sqrt(square(vec2[0]-vec1[0])+square(vec2[1]-vec1[1])+square(vec2[2]-vec1[2]));

	   
	}
	 
       }
     }
 

   ofstream Hneigh;
Hneigh.open("Hneigh.data");
Hneigh<<"S H N F FD      S   SD      OHO     ALEX<   MATTI   OO"<<endl;
  for(int i=0;i<nstruct/skip;i++)
  {      
    for(int j=0; j<nElt1; j++)
    {    //Hneigh<<j<<" "<<Onum[j+i*nElt1]<<" "<<Hbondneigh[j*OHsamples+i*nElt1*OHsamples]<<" "<<Hbondneighdist[j*OHsamples+i*nElt1*OHsamples]<<" ";
      for(int k=0; k<OHsamples;k++)
      {
	if(Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k] !=0 && Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k] != Hbondneighdist[j*OHsamples+i*nElt1*OHsamples])
	{
      Hneigh<<i<<" "<<j<<" "<<Onum[j+i*nElt1]<<" "<<Hbondneigh[j*OHsamples+i*nElt1*OHsamples]<<" "<<Hbondneighdist[j*OHsamples+i*nElt1*OHsamples]<<" "<<Hbondneigh[j*OHsamples+i*nElt1*OHsamples+k]<<" "<<Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k]<<" "<<anglesOHOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+(k)]*(180.0/3.14159)<<" " << 1.37+((-1.71)*(cos(anglesOHOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+(k)])) )<<" "<< anglesHOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+(k)]*(180.0/3.14159)<<" "<<distOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+(k)]<<endl;
	}
      //Hneigh<<j<<" "<<Onum[j+i*nElt1]<<" "<<Hbondneigh[j*OHsamples+i*nElt1*OHsamples]<<" "<<Hbondneighdist[j*OHsamples+i*nElt1*OHsamples]<<" "<<Hbondneigh[j*OHsamples+i*nElt1*OHsamples+1]<<" "<<Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+1]<<" "<<Hbondneigh[j*OHsamples+i*nElt1*OHsamples+2]<<" "<<Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+2]<<" "<<Hbondneigh[j*OHsamples+i*nElt1*OHsamples+3]<<" "<<Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+3]<<endl;
      }
    } 
  }  
Hneigh.close();
  
cout<<"Done!!!"<<endl;

/*
   ofstream Hneigh;
Hneigh.open("Hneigh.data");

  for(int i=0;i<nstruct/skip;i++)
  {      
    for(int j=0; j<nElt1; j++)
    {    
      Hneigh<<j<<" "<<Onum[j+i*nElt1]<<" "<<Hbondneigh[j*OHsamples+i*nElt1*OHsamples]<<" "<<Hbondneighdist[j*OHsamples+i*nElt1*OHsamples]<<" "<<Hbondneigh[j*OHsamples+i*nElt1*OHsamples+1]<<" "<<Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+1]<<" "<< Hbondxyz[j*5+i*nElt1*5+3]*(180.0/3.14159)<<" "<< angHoo[j+i*nElt1]*(180.0/3.14159)<<" "<<Hbondxyz[j*5+i*nElt1*5+4]<<endl;
      Hneigh<<j<<" "<<Onum[j+i*nElt1]<<" "<<Hbondneigh[j*OHsamples+i*nElt1*OHsamples]<<" "<<Hbondneighdist[j*OHsamples+i*nElt1*OHsamples]<<" "<<Hbondneigh[j*OHsamples+i*nElt1*OHsamples+1]<<" "<<Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+1]<<" "<<Hbondneigh[j*OHsamples+i*nElt1*OHsamples+2]<<" "<<Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+2]<<" "<<Hbondneigh[j*OHsamples+i*nElt1*OHsamples+3]<<" "<<Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+3]<<endl;
      
    } 
  }  
Hneigh.close();
*/

 ///Histograms here- yippie - I can derive the fit for Hbond criterium from here.
 
  int whichwater;
      if(HB_for=="int"){whichwater=1;}
      else if(HB_for=="bulk"){whichwater=1;}
      else if(HB_for=="all"){whichwater=2;} 
      
 if(HB_histograms=="yes")
 {
 
 cout<<"Histograms for angles and distance combinations "<<endl;
 
 float  OHint=0, OHmax=2;
 float  O_Hint=0, O_Hmax=bondist;
 float  OOint=0, OOmax=5.0;
 float  OHOint=0, OHOmax=180;
 float  HOOint=0, HOOmax=180;
 int num_binsOH=500;
 int num_binsO_H=500;
 int num_binsOO=500;
 int num_binsOHO=500;
 int num_binsHOO=500;
  OHint=OHmax/float(num_binsOH);
  O_Hint=O_Hmax/float(num_binsO_H);
  OOint=OOmax/float(num_binsOO);
  OHOint=OHOmax/float(num_binsOHO);
  HOOint=HOOmax/float(num_binsHOO);
    
  float *OHtics, *O_Htics, *OOtics, *OHOtics, *HOOtics;
  
  
  OHtics = (float *) malloc (sizeof(float)*num_binsOH*2);
  O_Htics = (float *) malloc (sizeof(float)*num_binsO_H*2);
  OOtics = (float *) malloc (sizeof(float)*num_binsOH*2);
  OHOtics = (float *) malloc (sizeof(float)*num_binsO_H*2);
  HOOtics = (float *) malloc (sizeof(float)*num_binsO_H*2);
  
  int *OH_O_H, *OH_OO, *OH_OHO, *OH_HOO;
  
  OH_O_H = (int *) malloc (sizeof(int)*num_binsOH*num_binsO_H);
  OH_OO = (int *) malloc (sizeof(int)*num_binsOH*num_binsOO);
  OH_OHO = (int *) malloc (sizeof(int)*num_binsOH*num_binsOHO);
  OH_HOO = (int *) malloc (sizeof(int)*num_binsOH*num_binsHOO);
  
  int *O_H_OO, *O_H_OHO, *O_H_HOO;
  
  O_H_OO = (int *) malloc (sizeof(int)*num_binsO_H*num_binsOO);
  O_H_OHO = (int *) malloc (sizeof(int)*num_binsO_H*num_binsOHO);
  O_H_HOO = (int *) malloc (sizeof(int)*num_binsO_H*num_binsHOO);
  
  int *OO_OHO, *OO_HOO;
  
  OO_OHO = (int *) malloc (sizeof(int)*num_binsOO*num_binsOHO);
  OO_HOO = (int *) malloc (sizeof(int)*num_binsOO*num_binsHOO);
  
  int *OHO_HOO;
  
  OHO_HOO = (int *) malloc (sizeof(int)*num_binsOHO*num_binsHOO);
  
  for(int i=0; i<num_binsOH; i++)
  {
   for(int j=0; j<num_binsO_H; j++)
   {
     OH_O_H[j+i*num_binsO_H]=0;
   }
   for(int j=0; j<num_binsOO; j++)
   {
     OH_OO[j+i*num_binsOO]=0;
   }
   for(int j=0; j<num_binsOHO; j++)
   {
     OH_OHO[j+i*num_binsOHO]=0;
   }
   for(int j=0; j<num_binsHOO; j++)
   {
     OH_HOO[j+i*num_binsHOO]=0;
   }
  }
  
  for(int i=0; i<num_binsO_H; i++)
  {
   for(int j=0; j<num_binsOO; j++)
   {
     O_H_OO[j+i*num_binsOO]=0;
   }
   for(int j=0; j<num_binsOHO; j++)
   {
     O_H_OHO[j+i*num_binsOHO]=0;
   }
   for(int j=0; j<num_binsHOO; j++)
   {
     O_H_HOO[j+i*num_binsHOO]=0;
   }
  }
  
  for(int i=0; i<num_binsOO; i++)
  {
   for(int j=0; j<num_binsOHO; j++)
   {
     OO_OHO[j+i*num_binsOHO]=0;
   }
   for(int j=0; j<num_binsHOO; j++)
   {
     OO_HOO[j+i*num_binsHOO]=0;
   }
  }
  
  for(int i=0; i<num_binsOHO; i++)
  {
   for(int j=0; j<num_binsHOO; j++)
   {
     OHO_HOO[j+i*num_binsHOO]=0;
   }
  }

  
for(int i=0; i<num_binsOH; i++)
{
 OHtics[i*2]=i*OHint;
 OHtics[i*2+1]=(i+1)*OHint;
}  
for(int i=0; i<num_binsO_H; i++)
{
 O_Htics[i*2]=i*O_Hint;
 O_Htics[i*2+1]=(i+1)*O_Hint;
} 
for(int i=0; i<num_binsOO; i++)
{
 OOtics[i*2]=i*OOint;
 OOtics[i*2+1]=(i+1)*OOint;
} 
for(int i=0; i<num_binsOHO; i++)
{
 OHOtics[i*2]=i*OHOint;
 OHOtics[i*2+1]=(i+1)*OHOint;
} 
for(int i=0; i<num_binsHOO; i++)
{
 HOOtics[i*2]=i*HOOint;
 HOOtics[i*2+1]=(i+1)*HOOint;
} 

/*
for(int i=0; i<num_binsOH; i++)
{
cout<< OHtics[i*2]<<" "<< OHtics[i*2+1]<<endl;
} 
for(int i=0; i<num_binsO_H; i++)
{
cout<< O_Htics[i*2]<<" "<< O_Htics[i*2+1]<<endl;
} 
*/
if(whichwater == 2)
{  
for(int i=0;i<nstruct/skip;i++)
  {
    for(int j=0; j<nElt1; j++)
    {
      //for(int k=0; k<OHsamples; k++)
      for(int k=0; k<5; k++) //only the four neighbouring water molecules
      {
	if(Onum[j+i*nElt1] >1 && Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k] > 0.0 && Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k] != Hbondneighdist[j*OHsamples+i*nElt1*OHsamples])
	{
         
	  for(int i1=0; i1<num_binsOH; i1++)
          {
	   if(Hbondneighdist[j*OHsamples+i*nElt1*OHsamples] > OHtics[i1*2] && Hbondneighdist[j*OHsamples+i*nElt1*OHsamples] < OHtics[i1*2+1])
	   {
	   
	    for(int i2=0; i2<num_binsO_H; i2++)
            {
             if(Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k] > O_Htics[i2*2] && Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k] < O_Htics[i2*2+1])
	     {
	      OH_O_H[i2+i1*num_binsO_H]+=1;
	      break;
	     }
	    }
	    for(int i2=0; i2<num_binsOO; i2++)
            {
             if(distOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k] > OOtics[i2*2] && distOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k] < OOtics[i2*2+1])
	     {
	      OH_OO[i2+i1*num_binsOO]+=1;
	      break;
	     }
	    }
	    for(int i2=0; i2<num_binsOHO; i2++)
            {
             if(anglesOHOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) > OHOtics[i2*2] && anglesOHOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) < OHOtics[i2*2+1])
	     {
	      OH_OHO[i2+i1*num_binsOHO]+=1;
	      break;
	     }
	    }
	    for(int i2=0; i2<num_binsHOO; i2++)
            {
             if(anglesHOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) > HOOtics[i2*2] && anglesHOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) < HOOtics[i2*2+1])
	     {
	      OH_HOO[i2+i1*num_binsHOO]+=1;
	      break;
	     }
	    }
	   break;
	   }//if for il ends
         }//i1 ends
         
         
         for(int i1=0; i1<num_binsO_H; i1++)
          {
	   if(Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k] > O_Htics[i1*2] && Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k] < O_Htics[i1*2+1])
	   {
	   
	    for(int i2=0; i2<num_binsOO; i2++)
            {
             if(distOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k] > OOtics[i2*2] && distOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k] < OOtics[i2*2+1])
	     {
	      O_H_OO[i2+i1*num_binsOO]+=1;
	      break;
	     }
	    }
	    for(int i2=0; i2<num_binsOHO; i2++)
            {
             if(anglesOHOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) > OHOtics[i2*2] && anglesOHOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) < OHOtics[i2*2+1])
	     {
	      O_H_OHO[i2+i1*num_binsOHO]+=1;
	      break;
	     }
	    }
	    for(int i2=0; i2<num_binsHOO; i2++)
            {
             if(anglesHOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) > HOOtics[i2*2] && anglesHOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) < HOOtics[i2*2+1])
	     {
	      O_H_HOO[i2+i1*num_binsHOO]+=1;
	      break;
	     }
	    }
	    break;
	   }//if for il ends
         }//i1 ends
         
         for(int i1=0; i1<num_binsOO; i1++)
          {
	   if(distOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k] > OOtics[i1*2] && distOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k] < OOtics[i1*2+1])
	   {
	   
	    for(int i2=0; i2<num_binsOHO; i2++)
            {
             if(anglesOHOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) > OHOtics[i2*2] && anglesOHOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) < OHOtics[i2*2+1])
	     {
	      OO_OHO[i2+i1*num_binsOHO]+=1;
	      break;
	     }
	    }
	    for(int i2=0; i2<num_binsHOO; i2++)
            {
             if(anglesHOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) > HOOtics[i2*2] && anglesHOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) < HOOtics[i2*2+1])
	     {
	      OO_HOO[i2+i1*num_binsHOO]+=1;
	      break;
	     }
	    }
	    break;
	   }//if for il ends
         }//i1 ends
         
         for(int i1=0; i1<num_binsOHO; i1++)
          {
	   if(anglesOHOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) > OHOtics[i1*2] && anglesOHOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) < OHOtics[i1*2+1])
	   {
	   
	    for(int i2=0; i2<num_binsHOO; i2++)
            {
             if(anglesHOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) > HOOtics[i2*2] && anglesHOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) < HOOtics[i2*2+1])
	     {
	      OHO_HOO[i2+i1*num_binsHOO]+=1;
	      break;
	     }
	    }
	   break;
	   }//if for il ends
         }//i1 ends
      }//if for k ends
     }//k ends
   }//j ends
  }//i ends
}  

else if(whichwater != 2)
{  
for(int i=0;i<nstruct/skip;i++)
  {
    for(int j=0; j<nElt1; j++)
    {
     if(exch1[j+i*nElt1] == whichwater)
     {
      //for(int k=0; k<OHsamples; k++)
       for(int k=0; k<5; k++) ///only the four neighbouring water molecules
      {
	if(Onum[j+i*nElt1] >1 && Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k] >0.0 && Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k] != Hbondneighdist[j*OHsamples+i*nElt1*OHsamples])
	{
         
	  for(int i1=0; i1<num_binsOH; i1++)
          {
	   if(Hbondneighdist[j*OHsamples+i*nElt1*OHsamples] > OHtics[i1*2] && Hbondneighdist[j*OHsamples+i*nElt1*OHsamples] < OHtics[i1*2+1])
	   {
	   
	    for(int i2=0; i2<num_binsO_H; i2++)
            {
             if(Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k] > O_Htics[i2*2] && Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k] < O_Htics[i2*2+1])
	     {
	      OH_O_H[i2+i1*num_binsO_H]+=1;
	      break;
	     }
	    }
	    for(int i2=0; i2<num_binsOO; i2++)
            {
             if(distOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k] > OOtics[i2*2] && distOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k] < OOtics[i2*2+1])
	     {
	      OH_OO[i2+i1*num_binsOO]+=1;
	      break;
	     }
	    }
	    for(int i2=0; i2<num_binsOHO; i2++)
            {
             if(anglesOHOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) > OHOtics[i2*2] && anglesOHOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) < OHOtics[i2*2+1])
	     {
	      OH_OHO[i2+i1*num_binsOHO]+=1;
	      break;
	     }
	    }
	    for(int i2=0; i2<num_binsHOO; i2++)
            {
             if(anglesHOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) > HOOtics[i2*2] && anglesHOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) < HOOtics[i2*2+1])
	     {
	      OH_HOO[i2+i1*num_binsHOO]+=1;
	      break;
	     }
	    }
	   break;
	   }//if for il ends
         }//i1 ends
         
         
         for(int i1=0; i1<num_binsO_H; i1++)
          {
	   if(Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k] > O_Htics[i1*2] && Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k] < O_Htics[i1*2+1])
	   {
	   
	    for(int i2=0; i2<num_binsOO; i2++)
            {
             if(distOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k] > OOtics[i2*2] && distOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k] < OOtics[i2*2+1])
	     {
	      O_H_OO[i2+i1*num_binsOO]+=1;
	      break;
	     }
	    }
	    for(int i2=0; i2<num_binsOHO; i2++)
            {
             if(anglesOHOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) > OHOtics[i2*2] && anglesOHOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) < OHOtics[i2*2+1])
	     {
	      O_H_OHO[i2+i1*num_binsOHO]+=1;
	      break;
	     }
	    }
	    for(int i2=0; i2<num_binsHOO; i2++)
            {
             if(anglesHOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) > HOOtics[i2*2] && anglesHOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) < HOOtics[i2*2+1])
	     {
	      O_H_HOO[i2+i1*num_binsHOO]+=1;
	      break;
	     }
	    }
	  break;
	   }//if for il ends
         }//i1 ends
         
         for(int i1=0; i1<num_binsOO; i1++)
          {
	   if(distOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k] > OOtics[i1*2] && distOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k] < OOtics[i1*2+1])
	   {
	   
	    for(int i2=0; i2<num_binsOHO; i2++)
            {
             if(anglesOHOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) > OHOtics[i2*2] && anglesOHOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) < OHOtics[i2*2+1])
	     {
	      OO_OHO[i2+i1*num_binsOHO]+=1;
	      break;
	     }
	    }
	    for(int i2=0; i2<num_binsHOO; i2++)
            {
             if(anglesHOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) > HOOtics[i2*2] && anglesHOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) < HOOtics[i2*2+1])
	     {
	      OO_HOO[i2+i1*num_binsHOO]+=1;
	      break;
	     }
	    }
	  break;
	   }//if for il ends
         }//i1 ends
         
         for(int i1=0; i1<num_binsOHO; i1++)
          {
	   if(anglesOHOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) > OHOtics[i1*2] && anglesOHOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) < OHOtics[i1*2+1])
	   {
	   
	    for(int i2=0; i2<num_binsHOO; i2++)
            {
             if(anglesHOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) > HOOtics[i2*2] && anglesHOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k]*(180.0/3.14159) < HOOtics[i2*2+1])
	     {
	      OHO_HOO[i2+i1*num_binsHOO]+=1;
	      break;
	     }
	    }
	    break;
	   }//if for il ends
         }//i1 ends
      }//if for k ends
     }//k ends
     }
   }//j ends
  }//i ends
}
  ofstream DOH_AOH, DOH_OO, DOH_OHO, DOH_HOO, AOH_OO, AOH_OHO, AOH_HOO, pOO_HOO, pOO_OHO, pOHO_HOO;
  
  DOH_AOH.open("DOH_AOH.data");
  DOH_OO.open("DOH_OO.data");
  DOH_OHO.open("DOH_OHO.data");
  DOH_HOO.open("DOH_HOO.data");
  AOH_OO.open("AOH_OO.data");
  AOH_OHO.open("AOH_OHO.data");
  AOH_HOO.open("AOH_HOO.data");
  pOO_OHO.open("OO_OHO.data");
  pOO_HOO.open("OO_HOO.data");
  pOHO_HOO.open("OHO_HOO.data");
  
  
  float normnum;
  
      if(HB_for=="int"){normnum = AVGintnum;}
      //if(HB_for=="int"){normnum = nElt1;}
      else if(HB_for=="bulk"){normnum = AVGbulknum;}
      else if(HB_for=="all"){normnum = AVGintnum+AVGbulknum;} 
      
      
  for(int i=0; i<num_binsOH; i++)
  {
   for(int j=0; j<num_binsO_H; j++)
   {
     DOH_AOH<<(OHtics[i*2]+OHtics[i*2+1])/2<<" "<<(O_Htics[j*2]+O_Htics[j*2+1])/2<<" "<<OH_O_H[j+i*num_binsO_H]/normnum<<endl;
   }DOH_AOH<<endl;
   for(int j=0; j<num_binsOO; j++)
   {
     DOH_OO<<(OHtics[i*2]+OHtics[i*2+1])/2<<" "<<(OOtics[j*2]+OOtics[j*2+1])/2<<" "<<OH_OO[j+i*num_binsOO]/normnum<<endl;
   } DOH_OO<<endl;
   for(int j=0; j<num_binsOHO; j++)
   {
     DOH_OHO<<(OHtics[i*2]+OHtics[i*2+1])/2<<" "<<(OHOtics[j*2]+OHOtics[j*2+1])/2<<" "<<OH_OHO[j+i*num_binsOHO]/normnum<<endl;
   }DOH_OHO<<endl;
   for(int j=0; j<num_binsHOO; j++)
   {
     DOH_HOO<<(OHtics[i*2]+OHtics[i*2+1])/2<<" "<<(HOOtics[j*2]+HOOtics[j*2+1])/2<<" "<<OH_HOO[j+i*num_binsHOO]/normnum<<endl;
   }DOH_HOO<<endl;
   
  }
  
  for(int i=0; i<num_binsO_H; i++)
  {
   for(int j=0; j<num_binsOO; j++)
   {
     AOH_OO<<(O_Htics[i*2]+O_Htics[i*2+1])/2<<" "<<(OOtics[j*2]+OOtics[j*2+1])/2<<" "<<O_H_OO[j+i*num_binsOO]/normnum<<endl;
   }AOH_OO<<endl;
   for(int j=0; j<num_binsOHO; j++)
   {
     AOH_OHO<<(O_Htics[i*2]+O_Htics[i*2+1])/2<<" "<<(OHOtics[j*2]+OHOtics[j*2+1])/2<<" "<<O_H_OHO[j+i*num_binsOHO]/normnum<<endl;
   }AOH_OHO<<endl;
   for(int j=0; j<num_binsHOO; j++)
   {
     AOH_HOO<<(O_Htics[i*2]+O_Htics[i*2+1])/2<<" "<<(HOOtics[j*2]+HOOtics[j*2+1])/2<<" "<<O_H_HOO[j+i*num_binsHOO]/normnum<<endl;
   }AOH_HOO<<endl;
  }
  
  for(int i=0; i<num_binsOO; i++)
  {
   for(int j=0; j<num_binsOHO; j++)
   {
     pOO_OHO<<(OOtics[i*2]+OOtics[i*2+1])/2<<" "<<(OHOtics[j*2]+OHOtics[j*2+1])/2<<" "<<OO_OHO[j+i*num_binsOHO]/normnum<<endl;
   }pOO_OHO<<endl;
   for(int j=0; j<num_binsHOO; j++)
   {
     pOO_HOO<<(OOtics[i*2]+OOtics[i*2+1])/2<<" "<<(HOOtics[j*2]+HOOtics[j*2+1])/2<<" "<<OO_HOO[j+i*num_binsHOO]/normnum<<endl;
   }pOO_HOO<<endl;
  }
  
  for(int i=0; i<num_binsOHO; i++)
  {
   for(int j=0; j<num_binsHOO; j++)
   {
     pOHO_HOO<<(OHOtics[i*2]+OHOtics[i*2+1])/2<<" "<<(HOOtics[j*2]+HOOtics[j*2+1])/2<<" "<<OHO_HOO[j+i*num_binsHOO]/normnum<<endl;
   }pOHO_HOO<<endl;
  }
  
  
  DOH_AOH.close();
  DOH_OO.close();
  DOH_OHO.close();
  DOH_HOO.close();
  AOH_OO.close();
  AOH_OHO.close();
  AOH_HOO.close();
  pOO_OHO.close();
  pOO_HOO.close();
  pOHO_HOO.close();

//Now, getting the Hbond per oxygen atom.

cout<<"Done !!!"<<endl;
}// histogram yes condition ends

ofstream Hbonds; //writes down the ids of H that are H bonded
Hbonds.open("Hbonds.data");
cout<<"Hydrogen bond identification based on the criterium set "<<endl;  
  
if(HB_criterium_set =="OOH-OO")  //matti
{
    for(int i=0;i<(nstruct/skip);i++)
     {Hbonds<<i<<endl;
      for(int j=0; j<nElt1; j++)
       {
	 for(int k=0; k<OHsamples; k++)
	 {
	   //actual H bond criterion
	   if(Onum[j+i*nElt1] >1 && Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k] !=0 && Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k] != Hbondneighdist[j*OHsamples+i*nElt1*OHsamples] && anglesHOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+(k)] <  Hbond_angle_dev && anglesHOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+(k)] !=  0.0&& distOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+(k)]<  max_O_O && distOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+(k)] != 0)
	   {
	    Hbondxyz[j*9+i*nElt1*9]=A[j*3+i*nElt1*3];
	    Hbondxyz[j*9+i*nElt1*9+1]=A[j*3+i*nElt1*3+1];
	    Hbondxyz[j*9+i*nElt1*9+2]=A[j*3+i*nElt1*3+2]; 
	    Hbondxyz[j*9+i*nElt1*9+3]= anglesHOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k];
	    Hbondxyz[j*9+i*nElt1*9+4]= distOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k];
	    Hbondxyz[j*9+i*nElt1*9+5]= Hbondneigh[j*OHsamples+i*nElt1*OHsamples+k];
	    Hbondxyz[j*9+i*nElt1*9+6]= Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k];
	    Hbondxyz[j*9+i*nElt1*9+7]= Hbondneighdist[j*OHsamples+i*nElt1*OHsamples];
	    Hbondxyz[j*9+i*nElt1*9+8]= anglesOHOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k];
	    Hbonds<<j<<" ";
	    break;
	     
	  }
	 }
       }
       Hbonds<<endl;
     }
}

Hbonds.close();

if(HB_criterium_set=="OHO-OH")   //alex
{
    for(int i=0;i<(nstruct/skip);i++)
     {
      for(int j=0; j<nElt1; j++)
       {
	 for(int k=0; k<OHsamples; k++)
	 {
	   if(Onum[j+i*nElt1] >1 && Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k] !=0 && Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k] != Hbondneighdist[j*OHsamples+i*nElt1*OHsamples] &&  Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k] != 0.0 &&  Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k] < 1.37+((-1.71)*(cos(anglesOHOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+(k)]))) )
	   {
	    Hbondxyz[j*9+i*nElt1*9]=A[j*3+i*nElt1*3];
	    Hbondxyz[j*9+i*nElt1*9+1]=A[j*3+i*nElt1*3+1];
	    Hbondxyz[j*9+i*nElt1*9+2]=A[j*3+i*nElt1*3+2];  
	    Hbondxyz[j*9+i*nElt1*9+3]= anglesOHOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k];
	    Hbondxyz[j*9+i*nElt1*9+4]= distOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k];
	    Hbondxyz[j*9+i*nElt1*9+5]= Hbondneigh[j*OHsamples+i*nElt1*OHsamples+k];
	    Hbondxyz[j*9+i*nElt1*9+6]= Hbondneighdist[j*OHsamples+i*nElt1*OHsamples+k];
	    Hbondxyz[j*9+i*nElt1*9+7]= Hbondneighdist[j*OHsamples+i*nElt1*OHsamples];
	    Hbondxyz[j*9+i*nElt1*9+8]= anglesHOOsamples[j*(OHsamples)+i*nElt1*(OHsamples)+k];
	    break;
	   }
	 }
       }
     }  
}
     




cout<<"Done !!!"<<endl;

cout<<"Mapping Hbonds to the donor acceptopr matrix for lifetime calculation"<<endl;
  
  
  ofstream Hexchplot, Hexchplotac;
 
 
 
 //copy coordinates for retrieval information
 int *initHBnum; // number of H bonds in each frame of the trajectory
 

//counting the number of Hbonds in each structure
    initHBnum = (int *) malloc (sizeof(int)*nstruct/skip);
    if(whichwater != 2)
    {
     for(int i=0;i<nstruct/skip;i++) //for each frame
     {
      initHBnum[i]=0;
      for(int j=0; j<nElt1; j++) // for each H atom
       {
	if(Hbondxyz[j*9+i*nElt1*9+4] !=0 && exch1[j+i*nElt1] == whichwater) 
        {
	   initHBnum[i]+=1;
	}
       }
     }
    }
    else if(whichwater ==2)
    {
     for(int i=0;i<nstruct/skip;i++)
     {
      initHBnum[i]=0;
      for(int j=0; j<nElt1; j++)
       { 
	if(Hbondxyz[j*9+i*nElt1*9+4] !=0) 
        {
	   initHBnum[i]+=1;
	}
       }
     }      
    }
    
//storing the coordinates of the Hbonds in nElt2*nElt2*nstruct/skip matrix

//hexch1 stores the acceptor - H matrix - for the H bond lifetime

//only H bond lifetime and no pair lifetimes!!!
if(HBlifetime == "yes")
   {

    float *Hexch, *Hexch1;
    Hexch = (float *) malloc (sizeof(float)*nElt2*nElt1*(nstruct/skip));
    Hexch1 = (float *) malloc (sizeof(float)*nElt2*nElt1*(nstruct/skip));
    
         for(int i=0;i<(nstruct/skip);i++)
         {
           for(int j=0; j<nElt2; j++)//acceptor oxygen
           {
	     for(int k=0; k<nElt1; k++)//shared hydrogen
             {
	       Hexch[k+j*nElt1+i*nElt2*nElt1]=0.0;
	       Hexch1[k+j*nElt1+i*nElt2*nElt1]=0.0;
	     }
	   }
	 }
cout <<"exch matrix"<<endl;	 

int **initHBcol1;// storing the location of the H bond in the donor/acceptor , acceptor-H matrix 

initHBcol1 = (int **) malloc (sizeof(int *)*nstruct/skip);
     for(int i=0; i<nstruct/skip; i++)// for each frame
     {
     initHBcol1[i] = (int *) malloc (sizeof(int)*initHBnum[i]);
     for(int j=0; j<initHBnum[i]; j++)// for number of H bonds this frame
     {
      initHBcol1[i][j]=0;
     }
     }
     
    if(whichwater != 2)
    {
    for(int i=0;i<(nstruct/skip);i++)
     {int chker=0;
      for(int j=0; j<nElt1; j++)
       {
	if( Hbondxyz[j*9+i*nElt1*9+4] !=0 && exch1[j+i*nElt1] == whichwater ) 
        {
	  Hexch1[j+Hbondneigh[j*OHsamples+i*nElt1*OHsamples+1]*nElt1+i*nElt2*nElt1]=1.0; //acceptor check
	  initHBcol1[i][chker] = (j+Hbondneigh[j*OHsamples+i*nElt1*OHsamples+1]*nElt1); chker+=1;
	}
	if( Hbondxyz[j*9+i*nElt1*9+4] !=0 ) 
        {
	  Hexch[j+Hbondneigh[j*OHsamples+i*nElt1*OHsamples+1]*nElt1+i*nElt2*nElt1]=1.0; //acceptor check
	}
	 
      }
     }
    }
    else if(whichwater == 2)
    {
     for(int i=0;i<(nstruct/skip);i++)
     {int chker=0;
      for(int j=0; j<nElt1; j++)
       {
	if(Hbondxyz[j*9+i*nElt1*9+4] !=0) 
        {
	  Hexch1[j+Hbondneigh[j*OHsamples+i*nElt1*OHsamples+1]*nElt1+i*nElt2*nElt1]=1.0;
	  Hexch[j+Hbondneigh[j*OHsamples+i*nElt1*OHsamples+1]*nElt1+i*nElt2*nElt1]=1.0;
	  initHBcol1[i][chker] = (j+Hbondneigh[j*OHsamples+i*nElt1*OHsamples+1]*nElt1); chker+=1;
	}
       }
     }  
    }
  /*
  lifetime
    Hexchplot.open("Hexch.data");

    for(int i=0;i<(nstruct/skip);i++)
         {
           for(int j=0; j<nElt2; j++)//donor
           {
	     for(int k=0; k<nElt2; k++)//acceptor
             {
	       Hexchplot << Hexch[k+j*nElt2+i*nElt2*nElt2]<<" ";
	     }
	     Hexchplot <<endl;
	   }
	   Hexchplot <<endl<<endl;
	 }
    
    Hexchplot.close(); */
cout<<" Done !!!"<<endl;     
  
  //HB lifetime check - timesaver
//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------  
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
 //########################################################################################################################################
//########################################################################################################################################
//BINNING AND AVERAGING THE HBAF
//#########################################################################################################################################
//######################################################################################################################################### 
	 
    for(int i=0; i<num_bins;i++)  //start bin
    {//int avgcounter=0; 
        printf("Cuda start: All lifetimes  %d \n",i);
     for(int j=0; j<SD_store; j=j+skips) //go through restart points
     // for(int j=0; j<1; j++)
     {      
       
       //cout<<initHBnum[j+i*SD_store]<<endl;
       float *Hexch_sized;//O-pair correlation
       float *Hexch_sized1;//H-bond correlation
       Hexch_sized = (float *) malloc (sizeof(float)*initHBnum[j+i*SD_store]*SD_store);
       Hexch_sized1 = (float *) malloc (sizeof(float)*initHBnum[j+i*SD_store]*SD_store);
         for(int ia=0;ia<SD_store;ia++)
         {
           for(int ja=0; ja<initHBnum[j+i*SD_store]; ja++)////only those H that are hydrogen bonded in the first frame in the restart bin
           { //if(i==0) cout << initHBcol[j]<<endl;
	     //cout<<ja<<" "<<initHBcol[j+i*SD_store][ja]<<endl;
	       Hexch_sized1[ja+ia*initHBnum[j+i*SD_store]] = Hexch1[initHBcol1[j+i*SD_store][ja]+(j+i*SD_store)*nElt2*nElt1+ia*nElt2*nElt1];
	       Hexch_sized[ja+ia*initHBnum[j+i*SD_store]]  = Hexch[initHBcol1[j+i*SD_store][ja]+(j+i*SD_store)*nElt2*nElt1+ia*nElt2*nElt1];
	   }
	 }
	 
     
     
 
      
       //SD=(float *) malloc (sizeof(float)*initHBnum[j+i*SD_store]*SD_store);
       SD1=(float *) malloc (sizeof(float)*initHBnum[j+i*SD_store]*SD_store);
        for(int ia=0;ia<SD_store;ia++)
        {
	 for(int ja=0; ja<initHBnum[j+i*SD_store]; ja++)
	 {
           //SD[ja+ia*initHBnum[j+i*SD_store]]=0.0;
	   SD1[ja+ia*initHBnum[j+i*SD_store]]=0.0;
	 }
	}
       
      
      dim3 dimBlocka(32,1,32);
      dim3 dimGrida((initHBnum[j+i*SD_store]+dimBlocka.x-1)/dimBlocka.x,1,(SD_store+dimBlocka.z-1)/dimBlocka.z);

       
      
      
      float *dev_A1;
            
      cudaMalloc((void **)&dev_A,sizeof(float)*initHBnum[j+i*SD_store]*SD_store);
      cudaMalloc((void **)&dev_A1,sizeof(float)*initHBnum[j+i*SD_store]*SD_store);
      cudaMalloc((void **)&dev_SD1,sizeof(float)*initHBnum[j+i*SD_store]*SD_store);
      
      cudaMemcpy(dev_A,Hexch_sized1,sizeof(float)*initHBnum[j+i*SD_store]*SD_store,cudaMemcpyHostToDevice);
      cudaMemcpy(dev_A1,Hexch_sized,sizeof(float)*initHBnum[j+i*SD_store]*SD_store,cudaMemcpyHostToDevice);
      cudaMemcpy(dev_SD1,SD1,sizeof(float)*initHBnum[j+i*SD_store]*SD_store,cudaMemcpyHostToDevice);
   
      HBAF_calc<<<dimGrida,dimBlocka>>>(dev_A,dev_A1,dev_SD1,SD_store,initHBnum[j+i*SD_store],i,j,origins,skips);
      cudaMemcpy(SD1,dev_SD1,sizeof(float)*initHBnum[j+i*SD_store]*SD_store,cudaMemcpyDeviceToHost);
      cudaFree(dev_A);
      cudaFree(dev_A1);
      
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
      SDreduce<<<dimGrids,dimBlocks>>>(dev_SD1,dev_SDavg1,SD_store,initHBnum[j+i*SD_store],fairy);
      cudaMemcpy(SDavg1,dev_SDavg1,sizeof(float)*SD_store,cudaMemcpyDeviceToHost);
       
      cudaFree(dev_SD1);
      cudaFree(dev_SDavg1);
      
           
      
      for(int ja =0; ja<SD_store;ja++)
      {
      //SDsum[ja]+=SDavg[ja];
      SDsum1[ja]+=SDavg1[ja];
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
      free(SD1);
      //avgcounter++;
    }
  
   // for(int ja =1; ja<SD_store-1;ja++)
   //   {
   //    FDsum[ja]/=avgcounter;
   //   }

    
   }
     //for(int ja =1; ja<SD_store;ja++)
     //{
      // SDsum[ja]/=SDsum[0];
      // SDsum1[ja]/=SDsum1[0];
     // } 
     // SDsum[0]=1.0;
     // SDsum1[0]=1.0;
      //Hexchplot.open("ct-pair.data");
      Hexchplotac.open("ct.data");
     for(int ja =0; ja<SD_store;ja++)
     {
       //SDsum[ja]/=float(num_bins*restarts);
       SDsum1[ja]/=float(num_bins*restarts);
       //Hexchplot <<ja*timestep/1000<<" "<<SDsum[ja]/SDsum[0]<<endl;
       Hexchplotac <<ja*timestep/1000<<" "<<SDsum1[ja]/SDsum1[0]<<endl;
       //Hexchplot <<ja*timestep/1000<<" "<<SDsum[ja]<<endl;
      } 
      //Hexchplot.close();
      Hexchplotac.close();
      
      //averaging for a proper plot
    
      /*
      for(int ja =0; ja<SD_store;ja++)
     {
       SDsum[ja]-=SDsum[SD_store-1];
       SDsum1[ja]-=SDsum1[SD_store-1];
      }*/
      
/*   
     cout << "Gnuplot fitting begins"<< endl;
  FILE *pipe_gp = popen("gnuplot", "w");
   fputs("set fit logfile 'fit-pair.data'\n",pipe_gp);
   fputs("set fit quiet\n",pipe_gp);
  fputs("f(x)=A*exp(-x/i)+B*exp(-x/j)+C*exp(-x/k)\n", pipe_gp);
  fputs("A               = 0.339959\n", pipe_gp);
  fputs("i               = 2.14969\n", pipe_gp);
  fputs("B               = 0.303315\n", pipe_gp);
  fputs("j               = 0.0457113\n", pipe_gp);
  fputs("C               = 0.259671\n", pipe_gp);
  fputs("k               = 7.19704\n", pipe_gp);
  fputs("fit f(x) 'ct-pair.data' u 1:2 via A, i, B, j, C, k\n", pipe_gp);
  fputs("exit\n", pipe_gp);
  pclose(pipe_gp);
  
  FILE *pipe_gpa = popen("gnuplot", "w");
fputs("set fit logfile 'fit.data'\n",pipe_gpa);
fputs("set fit quiet\n",pipe_gpa);
  fputs("f(x)=A*exp(-x/i)+B*exp(-x/j)+C*exp(-x/k)\n", pipe_gpa);
  fputs("A               = 0.339959\n", pipe_gpa);
  fputs("i               = 2.14969\n", pipe_gpa);
  fputs("B               = 0.303315\n", pipe_gpa);
  fputs("j               = 0.0457113\n", pipe_gpa);
  fputs("C               = 0.259671\n", pipe_gpa);
  fputs("k               = 7.19704\n", pipe_gpa);
  fputs("fit f(x) 'ct.data' u 1:2 via A, i, B, j, C, k\n", pipe_gpa);
  fputs("exit\n", pipe_gpa);
  pclose(pipe_gpa);
  
  FILE *pipe_gp2 = popen("gnuplot", "w");
fputs("set fit logfile 'fit2.data'\n",pipe_gp2);
fputs("set fit quiet\n",pipe_gp2);
  fputs("f(x)=A*exp(-x/i)+B*exp(-x/j)\n", pipe_gp2);
  fputs("A               = 0.326531\n", pipe_gp2);
  fputs("i               = 0.514234\n", pipe_gp2);
  fputs("B               = 0.592641\n", pipe_gp2);
  fputs("j               = 4.81101\n", pipe_gp2); 
  fputs("fit f(x) 'ct.data' u 1:2 via A, i, B, j\n", pipe_gp2);
  fputs("exit\n", pipe_gp2);
  pclose(pipe_gp2);
     
  FILE *pipe_gp2a = popen("gnuplot", "w");
fputs("set fit logfile 'fit2-pair.data'\n",pipe_gp2a);
fputs("set fit quiet\n",pipe_gp2a);
  fputs("f(x)=A*exp(-x/i)+B*exp(-x/j)\n", pipe_gp2a);
  fputs("A               = 0.326531\n", pipe_gp2a);
  fputs("i               = 0.514234\n", pipe_gp2a);
  fputs("B               = 0.592641\n", pipe_gp2a);
  fputs("j               = 4.81101\n", pipe_gp2a); 
  fputs("fit f(x) 'ct-pair.data' u 1:2 via A, i, B, j\n", pipe_gp2a);
  fputs("exit\n", pipe_gp2a);
  pclose(pipe_gp2a);
  
cout<<"Gnuplot fits done"<<endl;  
*/
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
       int_SDsum1[ja]= (((ja*timestep*0.001) - ((ja-1)*timestep*0.001)) * ((SDsum1[ja] + SDsum1[ja-1])/2.0));
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
     Hexchplot.open("Hbond-lifetime.data");
     Hexchplot <<"Lifetime from the integral of c(t) is : "<<lifetime1<<" ps"<<endl;
     Hexchplot.close();
     
 
 
 
    float *Hbondnum;
    
    Hbondnum =(float*) malloc (sizeof(float)*(nstruct/skip));
    for(int i=0; i<(nstruct/skip); i++)
    {
      Hbondnum[i]=0;
    }
    for(int i=0; i<(nstruct/skip); i++)
    {
      for(int j=0; j<nElt2; j++)//acceptor
      {
	for(int k=0; k<nElt1; k++)//shared hydrogen
	{
	  Hbondnum[i]+=Hexch1[k+j*nElt1+i*nElt2*nElt1];
	}
      }
      //cout << Hbondnum[i]<<endl;
      //Hbondnum[i]/=(0.5*nElt2*(nElt2-1));
    }
    float avg_Hbondnum=0.0;
    
        Hexchplot.open("Hbondnum.data");

         for(int i=0;i<(nstruct/skip);i++)
         {
	   avg_Hbondnum+=Hbondnum[i];
	   Hexchplot <<i<<" "<<Hbondnum[i]<<endl;
	 }
	 avg_Hbondnum /= (nstruct/skip);
        Hexchplot.close();
	Hexchplot.open("Avg_Hbondnum.data");
	Hexchplot <<"Average number of Hbonds per structure = "<<avg_Hbondnum<<endl;
        Hexchplot.close();
 
  free(Hexch);
  free(Hexch1);
  free(Hbondnum);
  free(initHBcol1);
}  
  free(A);
  free(B);
  free(Onum);
  free(Hbondneigh);
  free(Hbondneighdist);
 
  /////////////////Histogramming part - trivial/////////////////////////////////////////////////////////
  
cout<<"Hbonds and angles: " << nElt1 <<" "<< Elt1 << " atoms are there in each structure"<<endl;

int Dirn,split;
float min,max;

if(set_max_z == "yes")
{
  minz = set_minz;
  maxz = set_maxz;
  Dirn = 2;
  min=minz;
  max=maxz;
  split=zsplit;
  
}

else if(set_max_y == "yes")
{
  minz = set_miny;
  maxz = set_maxy;
  Dirn=1;
  min=miny;
  max=maxy;
  split=ysplit;
}
else if(set_max_x == "yes")
{
  minz = set_minx;
  maxz = set_maxx;
  Dirn=0;
  min=minx;
  max=maxx;
  split=xsplit;
}

printf("%f %f\n",min,max);

    zrange=max-min;

    zint = zrange/split;

    ztick=(float*) malloc (sizeof(float)*split*2);
    
    for(int i=0;i<split;i++)
    {
      ztick[i*2]=min+(i*zint);
      ztick[i*2+1]=min+((i+1)*zint);
    }
    

FILE *xyztick = fopen("xyztick.dat","wt");
    for(int i=0;i<split;i++)
    {
      fprintf(xyztick,"%f %f \n",ztick[i*2],ztick[i*2+1]);
    }
fclose(xyztick);

   float *Hbond_density, *dev_Hbond_density;
   float *angle_density, *dev_angle_density;
   float *distance_density, *dev_distance_density;
   float *OH_density, *dev_OH_density;
   float *O_H_density, *dev_O_H_density;
   float *ang_density, *dev_ang_density;

    Hbond_density=(float *) malloc (sizeof(float)*split);
    angle_density=(float *) malloc (sizeof(float)*split);
    distance_density=(float *) malloc (sizeof(float)*split);
    OH_density=(float *) malloc (sizeof(float)*split);
    O_H_density=(float *) malloc (sizeof(float)*split);
    ang_density=(float *) malloc (sizeof(float)*split);

        for(int k=0;k<split;k++)
        {
          Hbond_density[k]=0.0;
          angle_density[k]=0.0;
          distance_density[k]=0.0;
	  OH_density[k]=0.0;
          ang_density[k]=0.0;
          O_H_density[k]=0.0;
        }


    printf("Start of cuda calculation\n");
    
    cudaMalloc((void **)&dev_A,sizeof(float)*(nstruct/skip)*nElt1*9);
    cudaMalloc((void **)&dev_Hbond_density,sizeof(float)*split);
    cudaMalloc((void **)&dev_angle_density,sizeof(float)*split);
    cudaMalloc((void **)&dev_distance_density,sizeof(float)*split);
    cudaMalloc((void **)&dev_OH_density,sizeof(float)*split);
    cudaMalloc((void **)&dev_O_H_density,sizeof(float)*split);
    cudaMalloc((void **)&dev_ang_density,sizeof(float)*split);
    cudaMalloc((void **)&dev_ztick, sizeof(float)*split*2);
    
    cudaMemcpy(dev_A,Hbondxyz,sizeof(float)*(nstruct/skip)*nElt1*9,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_Hbond_density,Hbond_density,sizeof(float)*split,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_angle_density,angle_density,sizeof(float)*split,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_distance_density,distance_density,sizeof(float)*split,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_OH_density,OH_density,sizeof(float)*split,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_ang_density,ang_density,sizeof(float)*split,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_O_H_density,O_H_density,sizeof(float)*split,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_ztick,ztick,sizeof(float)*split*2,cudaMemcpyHostToDevice);
    
    cout << "Memory for storing structure data: "<<((nstruct/skip)*nElt1*9)*4/float(1000000000)<< " Gbs" << endl;
    cout << "Memory for local storage: "<<((split*6)+(split*2))*4/float(1000000000)<< " Gbs" << endl;   
    
    dim3 dimBlockz(32,32,1);
    //dim3 dimBlock(1,1,1);
    dim3 dimGridz(((nstruct/skip)+dimBlockz.x-1)/dimBlockz.x,(nElt1+dimBlockz.y-1)/dimBlockz.y,1);
    //dim3 dimGrid(nstruct,nElt1,1);
    Hbond_calc<<<dimGridz,dimBlockz>>>(dev_A,dev_Hbond_density,dev_angle_density,dev_distance_density,dev_OH_density,dev_O_H_density,dev_ang_density,(nstruct/skip),nElt1,split,dev_ztick,Hbond_angle_dev,max_O_O,Dirn);
    cudaMemcpy(Hbond_density,dev_Hbond_density,sizeof(float)*split,cudaMemcpyDeviceToHost);
    cudaMemcpy(angle_density,dev_angle_density,sizeof(float)*split,cudaMemcpyDeviceToHost);
    cudaMemcpy(distance_density,dev_distance_density,sizeof(float)*split,cudaMemcpyDeviceToHost);
    cudaMemcpy(OH_density,dev_OH_density,sizeof(float)*split,cudaMemcpyDeviceToHost);
    cudaMemcpy(O_H_density,dev_O_H_density,sizeof(float)*split,cudaMemcpyDeviceToHost);
    cudaMemcpy(ang_density,dev_ang_density,sizeof(float)*split,cudaMemcpyDeviceToHost);
  
    float cellvolume;
    
    if(cell_type=="orthorhombic")
    {
    if(set_max_z == "yes") cellvolume = xvec*yvec*zint;
      if(set_max_x == "yes") cellvolume = yvec*zvec*zint;
      if(set_max_y == "yes") cellvolume = xvec*zvec*zint;
    }
    else if(cell_type=="monoclinic")
    {
    if(set_max_z == "yes") cellvolume= xvec*xvec*zint*0.866025403; //*sin(60) for rhombohedral 
    if(set_max_x == "yes") cellvolume= yvec*zvec*zint*0.866025403;
    if(set_max_y == "yes") cellvolume= xvec*zvec*zint*0.866025403;
      
    }
    
FILE *Zdistfile=fopen("Hbonddist.data","wt");

int zstart = 0;
int zend = split;

for(int m=zstart;m<zend;m++)
{
fprintf(Zdistfile,"%f %f %f %f %f %f %f \n ",ztick[m*2+1], Hbond_density[m]/(cellvolume*(nstruct/skip)),
        (angle_density[m]/Hbond_density[m])*(180.0/3.14159),distance_density[m]/Hbond_density[m], OH_density[m]/Hbond_density[m], O_H_density[m]/Hbond_density[m], (ang_density[m]/Hbond_density[m])*(180.0/3.14159));
}
fclose(Zdistfile);

cudaFree(dev_Hbond_density);
cudaFree(dev_A);
cudaFree(dev_ztick);
cudaFree(dev_distance_density);
cudaFree(dev_angle_density);
cudaFree(dev_OH_density);
cudaFree(dev_ang_density);
cudaFree(dev_O_H_density);

free(Hbondxyz);
free(Hbond_density);
free(angle_density);
free(distance_density);
free(OH_density);
free(ang_density);
free(O_H_density);
free(ztick);
  
  
  
  
}
