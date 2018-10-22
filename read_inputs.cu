//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"


void read_inputs(void)
{
 int lnElt1,lnElt2;
 nstruct=0;
 lnElt2=0;
 lnElt1=0;
 
// Read the number of structures and atoms per structure in the trajectory
 //command = "grep -c -E 'Atoms|generated' input.xyz > nstruct";
 string cmd("grep -c -E 'Atoms|generated' ");
 cmd += file_path;
 cmd += " > nstruct";
 system (cmd.c_str());
 ifstream nstructure;
 nstructure.open("nstruct");
 nstructure >> nstruct;
 nstructure.close();
 system("rm nstruct");
 ifstream traj;
 if(max_frames!=0){nstruct=max_frames;}
 nstruct=nstruct-start_frame;
 //if(unwrap == "yes") {traj.open("input.xyz_unwrapped",ios::in);}
 //else {traj.open("input.xyz",ios::in);}
 cout <<"using file in "<< file_path<<endl;
 traj.open(file_path.c_str());
//count number of atoms in each frame
 traj>>natoms;
 traj.seekg(0);
 cout << natoms <<" atoms found in "<< nstruct <<" frames"<<endl; 

//Count the number of atoms per chosen element from the trajectory's first frame
 float *maxmin;
 maxmin= (float *) malloc (sizeof(float)*natoms*3);
 count_metal=0;
  for(int i = 0; i< 1; i++)
  {
    do
    traj.get(c);
    while ( c != '\n');
    do
    traj.get(c);
    while ( c != '\n');
    cout<<"counting atoms of chosen element"<<endl;
    for(int j=0; j<natoms; j++)
    {
      traj >> felt >> maxmin[j*3+i*natoms*3] >> maxmin[j*3+i*natoms*3+1]
           >> maxmin[j*3+i*natoms*3+2];
      if(felt == Elt1)
      {
        lnElt1+=1;
      }
      else if(felt == Elt2)
      {
        lnElt2+=1;
      }
      else if(felt == rdf_metal_exclude)
      {
	count_metal+=1;
      }
      else if (Elt1 == "water")
      {
	if(felt =="O" || felt =="H")
	{
	  lnElt1+=1;
	}
      }
      if(felt == metal_species)
      {
	count_metal+=1;
      }
      do
      traj.get(c);
      while ( c != '\n');
    }
  }

  traj.seekg(0);
  if(Elt1==Elt2)
  {
    lnElt2=lnElt1;
    cout <<Elt1<<" = "<<lnElt1<<endl;
  }
  else
  {
    cout <<Elt1 <<" = "<< lnElt1<<" ; "<<Elt2<<" = "<<lnElt2<<endl;    
  }
  if(Elt1=="All" && Elt2=="All") {lnElt1=natoms; lnElt2=natoms;}
  if(type=="rdf") count_metal/=2;
//Now to map the volume spanned by the selected atom species,
//the minimum and maximum values of x,y,z coordinates are
//computed and stored - useful for the density/velocity plots.
int stride=0;
    minx = minimum(maxmin,natoms,stride);
    miny = minimum(maxmin,natoms,stride+1);
    minz = minimum(maxmin,natoms,stride+2);
    maxx = maximum(maxmin,natoms,stride);
    maxy = maximum(maxmin,natoms,stride+1);
    maxz = maximum(maxmin,natoms,stride+2);
    free(maxmin);

    
    
//Allocate matrices for the element coordinates as well as the 
//distance and histogram components
//Matrix A contains positions of species 1 and in case of RDF, B matrix
//holds the positions of the pair element.
//In case of same species distribution function, A and B are identical.

  if(average=="yes")
  {
   Aavg= (float *) malloc (sizeof(float)*lnElt1*3);
   for(int i=0; i<lnElt1*3;i++)
   {
    Aavg[i]=0.0;
   }
  }

      //experimental ///mainly for the diffusion density of int/bulk water
if(scope =="bulk-interface")  
  {
    float *first_frame, *metal_frame;
    first_frame=(float*) malloc(sizeof(float)*natoms*3);
    metal_frame=(float*) malloc(sizeof(float)*count_metal*3);
    int *species_order;
    species_order=(int *) malloc (sizeof(int)*natoms);
    for(int i=0;i<natoms;i++)
    {
      species_order[i]=0;
    }
    //Read coordinates of respective elements into the matrices
printf("Reading the positions of first frame in to data structures\n");
float Dum1,Dum2,Dum3;
//int *globalindexA, *globalindexB, *globalindexM;
//globalindexA = (int *) malloc (sizeof(int) * lnElt1);
//globalindexB = (int *) malloc (sizeof(int) * lnElt2);
//globalindexm = (int *) malloc (sizeof(int) * count_metal);
 for(int i=0;i<(start_frame*(natoms+2));i++)
 {
    getline(traj,Dummy);
 }
 for(int i=0; i< 1;i=i+1)
  {
    //int counta=0, countb=0, 
    int countc=0;
    getline(traj,Dummy);
    getline(traj,Dummy);
    for(int j=0; j< natoms; j++)
    {
      traj >> felt >> Dum1 >> Dum2 >> Dum3;
      first_frame[j*3]=Dum1;
      first_frame[j*3+1]=Dum2;
      first_frame[j*3+2]=Dum3;
      if(felt==Elt1){species_order[j]=1;}
      else if(felt==Elt2){species_order[j]=2;}
      /*if(felt==Elt1)
      {
        A[counta*3+i*lnElt1*3]= Dum1;
        A[counta*3+i*lnElt1*3+1]= Dum2;
        A[counta*3+i*lnElt1*3+2]= Dum3;
	globalindexA[counta]=j;
	counta=counta+1;
      }
      else if(felt==Elt2)
      {
        B[countb*3+i*lnElt2*3] = Dum1;
        B[countb*3+i*lnElt2*3+1] = Dum2;
        B[countb*3+i*lnElt2*3+2] = Dum3;
	globalindexB[countb]=j;
        countb=countb+1;
      }*/
      if(felt==metal_species)
      {
        metal_frame[countc*3] = Dum1;
        metal_frame[countc*3+1] = Dum2;
        metal_frame[countc*3+2] = Dum3;
	//globalindexM[countc]=j;
        countc=countc+1;
      }
    }
    getline(traj,Dummy);
    for (int m = 0 ; m < skip-1; m++)
    {
      for(int j=0;j<natoms;j++)
      {
	getline(traj,Dummy);
      }
      getline(traj,Dummy);
      getline(traj,Dummy);
    }
  }
  traj.seekg(0);
 // read the lattice information from the last line of the trajectory
 string cmd("tail -n 1 ");
 cmd += file_path;
 cmd += " > cellsize";
 system (cmd.c_str());
 ifstream cellsize;
 cellsize.open("cellsize",ios::in);
 yxvec =0; 
 if(cell_type == "orthorhombic") cellsize>> xvec >>yvec>> zvec;
 else if(cell_type == "monoclinic") cellsize >>xvec>>yvec>>zvec>>yxvec;
 cellsize.close();
 system("rm cellsize");
 cout <<"File read complete"<<endl;
 cout <<"Lattice length (X Y Z) = "<< xvec << " "<< yvec << " "<< zvec << " " << yxvec <<endl;
 
    
    nA_int=0;
    nB_int=0;
    nA_bulk=0;
    nB_bulk=0;
    int *local_distmat;
    local_distmat = (int *) malloc (sizeof(int)*natoms);
      for(int j=0; j<natoms; j++)
      {
	local_distmat[j]=0;
      }

//int *localindexA, *localindexB, *localindexM;
//localindexA = (int *) malloc (sizeof(int) * lnElt1);
//localindexB = (int *) malloc (sizeof(int) * lnElt2);
//localindexm = (int *) malloc (sizeof(int) * count_metal); 

     for(int i=0; i<count_metal; i++)
    {
      for(int j=0; j<natoms; j++)
      {
	float chkx,chky,chkz;
	chkx=fabs(metal_frame[i*3]-first_frame[j*3]);
        chky=fabs(metal_frame[i*3+1]-first_frame[j*3+1]);
        chkz=fabs(metal_frame[i*3+2]-first_frame[j*3+2]);
	if(cell_type=="monoclinic")
	{
        if(chky > yvec/2) {chky = chky - yvec; chkx = chkx - yxvec;}
        if(chkx > xvec/2) {chkx = chkx - xvec;}
        if(chkz > zvec/2) {chkz = chkz - zvec;}
	if (sqrt((chkx*chkx) + (chky*chky) + (chkz*chkz)) < bondist ){local_distmat[j]+=1;}
	}
	else if(cell_type == "orthorhombic")
	{
        if(chky > yvec/2) {chky = chky - yvec;}
        if(chkx > xvec/2) {chkx = chkx - xvec;}
        if(chkz > zvec/2) {chkz = chkz - zvec;}
	if (sqrt((chkx*chkx) + (chky*chky) + (chkz*chkz)) < bondist ){local_distmat[j]+=1;}
	}
      }
    }
      for(int j=0; j<natoms; j++)
      {
	if(local_distmat[j] >=1 && species_order[j]==1){nA_int++;}
	else if(local_distmat[j] >=1 && species_order[j]==2){nB_int++;}
	if(local_distmat[j] <1 && species_order[j]==1){nA_bulk++;}
	else if(local_distmat[j] <1 && species_order[j]==2){nB_bulk++;}
	
      }


A_int=(float*)malloc(sizeof(float)*nA_int*(nstruct/skip)*3);
B_int=(float*)malloc(sizeof(float)*nB_int*(nstruct/skip)*3);
A_bulk=(float*)malloc(sizeof(float)*nA_bulk*(nstruct/skip)*3);
B_bulk=(float*)malloc(sizeof(float)*nB_bulk*(nstruct/skip)*3);

printf("Reading the positions in to data structures\n");
cout << "for every frame "<< skip-1 <<" frames will be skipped"<<endl;
cout << nstruct/skip <<" frames will be used for the analysis"<<endl;

 for(int i=0;i<(start_frame*(natoms+2));i++)
 {
    getline(traj,Dummy);
 }

 for(int i=0; i< nstruct/skip;i=i+1)
  {
    int counta=0, countb=0, counta1=0, countb1=0;
    getline(traj,Dummy);
    getline(traj,Dummy);
    for(int j=0; j< natoms; j++)
    {
      traj >> felt >> Dum1 >> Dum2 >> Dum3;
      if(felt==Elt1 && local_distmat[j] >=1)
      {
        A_int[counta*3+i*nA_int*3] = Dum1;
        A_int[counta*3+i*nA_int*3+1] = Dum2;
        A_int[counta*3+i*nA_int*3+2]= Dum3;
	counta=counta+1;
      }
      else if(felt==Elt2 && local_distmat[j] >=1)
      {
        B_int[countb*3+i*nB_int*3] = Dum1;
        B_int[countb*3+i*nB_int*3+1] = Dum2;
        B_int[countb*3+i*nB_int*3+2] = Dum3;
        countb=countb+1;
      }
      else if(felt==Elt1 && local_distmat[j] <1)
      {
        A_bulk[counta1*3+i*nA_bulk*3] = Dum1;
        A_bulk[counta1*3+i*nA_bulk*3+1] = Dum2;
        A_bulk[counta1*3+i*nA_bulk*3+2]= Dum3;
	counta1=counta1+1;
      }
      else if(felt==Elt2 && local_distmat[j] <1)
      {
        B_bulk[countb1*3+i*nB_bulk*3] = Dum1;
        B_bulk[countb1*3+i*nB_bulk*3+1] = Dum2;
        B_bulk[countb1*3+i*nB_bulk*3+2] = Dum3;
        countb1=countb1+1;
      }
    }
    getline(traj,Dummy); 
    for (int m = 0 ; m < skip-1; m++)
    {
      for(int j=0;j<natoms;j++)
      {
	getline(traj,Dummy);
      }
      getline(traj,Dummy);
      getline(traj,Dummy);
    }
  }
   traj.close();  
   

   
}



     //for everything else...
else if (scope =="bulk")
{
    //
  A= (float *) malloc (sizeof(float)*lnElt1*(nstruct/skip)*3);

  if (Elt2 != "")
  {
  B= (float *) malloc (sizeof(float)*lnElt2*(nstruct/skip)*3);
  }
  if (metal_species != "")
  {
    METAL= (float *) malloc (sizeof(float)*count_metal*(nstruct/skip)*3);
  }
  if (unwrap == "yes")
  {
      Aint= (float *) malloc (sizeof(float)*lnElt1*(nstruct/skip)*3);
  }
  
  if(Elt1=="All" && Elt2=="All")
  {
  A_all= (float *) malloc (sizeof(float)*lnElt1*(nstruct/skip)*4);
  B_all= (float *) malloc (sizeof(float)*lnElt2*(nstruct/skip)*4);
  }
//Read coordinates of respective elements into the matrices

 //cout << (nstruct/skip)<<endl;
  lattice = (float *) malloc (sizeof(float)*6*(nstruct/skip));
  for(int i =0; i<(nstruct/skip); i++)
  {
      for(int j=0; j<6;j++)
      {
      lattice[j+i*6]=0.0;
      }
  }


printf("Reading the positions in to data structures\n");
cout << "for every frame "<< skip-1 <<" frames will be skipped"<<endl;
cout << nstruct/skip <<" frames will be used for the analysis"<<endl;
float Dum1,Dum2,Dum3;

 for(int i=0;i<(start_frame*(natoms+2));i++)
 {
    getline(traj,Dummy);
 }
float Dum4;
 for(int i=0; i<(nstruct/skip);i=i+1)
  {// cout<<i<<endl; 
    int counta=0, countb=0, countc=0;
    getline(traj,Dummy);
    //getline(traj,Dummy);
    traj >> lattice[i*6] >> lattice[1+i*6] >> lattice[2+i*6]>> lattice[3+i*6]>>lattice[4+i*6]>>lattice[5+i*6];
    //cout <<" "<< lattice[i*6] <<" "<< lattice[1+i*6] <<" "<< lattice[2+i*6]<<endl;
    //getline(traj,Dummy);
    for(int j=0; j< natoms; j++)
    {
      
      if(ext_charge=="yes"){traj >> felt >> Dum1 >> Dum2 >> Dum3 >> Dum4;}
      else {traj >> felt >> Dum1 >> Dum2 >> Dum3;}
      if(felt==Elt1)
      {
        A[counta*3+i*lnElt1*3] = Dum1;
        A[counta*3+i*lnElt1*3+1] = Dum2;
        A[counta*3+i*lnElt1*3+2]= Dum3;
	counta=counta+1;
      }
      else if(felt==Elt2)
      {
        B[countb*3+i*lnElt2*3] = Dum1;
        B[countb*3+i*lnElt2*3+1] = Dum2;
        B[countb*3+i*lnElt2*3+2] = Dum3;
        countb=countb+1;
      }
      
      else if(Elt1=="water")
      {
	if(felt=="O" || felt=="H")
	{
	A[counta*3+i*lnElt1*3] = Dum1;
        A[counta*3+i*lnElt1*3+1] = Dum2;
        A[counta*3+i*lnElt1*3+2]= Dum3;
	counta=counta+1;
	}
      }
      if(felt==metal_species)
      {
        METAL[countc*3+i*count_metal*3] = Dum1;
        METAL[countc*3+i*count_metal*3+1] = Dum2;
        METAL[countc*3+i*count_metal*3+2] = Dum3;
        countc=countc+1;
      }
      if(Elt1=="All" && Elt2=="All")
      {
	A_all[counta*4+i*lnElt1*4] = Dum1;
        A_all[counta*4+i*lnElt1*4+1] = Dum2;
        A_all[counta*4+i*lnElt1*4+2]= Dum3;
	if(ext_charge=="yes"){A_all[counta*4+i*lnElt1*4+3] = Dum4;}
	else{
	if(felt=="O") {A_all[counta*4+i*lnElt1*4+3]= -2.0;}
	else if(felt=="H") {A_all[counta*4+i*lnElt1*4+3]= 1.0;}
	else if(felt=="Cu") {A_all[counta*4+i*lnElt1*4+3]= 0.0;}
	else if(felt=="Zn") {A_all[counta*4+i*lnElt1*4+3]= 2.0;}
	}
	counta=counta+1;
      }
      
    }
    getline(traj,Dummy);
    for (int m = 0 ; m < skip-1; m++)
    {
      for(int j=0;j<natoms;j++)
      {
	getline(traj,Dummy);
      }
      getline(traj,Dummy);
      getline(traj,Dummy);
    }
  }
  if(Elt1=="All" && Elt2=="All") {B_all = A_all;}
   traj.close();  
   //cout <<"HI"<<endl;
// read the lattice information from the last line of the trajectory
 string cmd("tail -n 1 ");
 cmd += file_path;
 cmd += " > cellsize";
 system (cmd.c_str());
 ifstream cellsize;
 cellsize.open("cellsize",ios::in);
 yxvec =0; 
 if(cell_type == "orthorhombic") cellsize>> xvec >>yvec>> zvec;
 else if(cell_type == "monoclinic") cellsize >>xvec>>yvec>>zvec>>yxvec;
 cellsize.close();
 system("rm cellsize");
 cout <<"File read complete"<<endl;
 cout <<"Lattice length (X Y Z) = "<< xvec << " "<< yvec << " "<< zvec << " " << yxvec <<endl;  
}
//velocity information is also read and stored incase of velocity plot
  if(veloc=="yes")
  {
    float Dum1,Dum2,Dum3;
  printf("Reading the velocities in to data structures\n");
  if(Elt1=="all")
  {
  VEL= (float *) malloc (sizeof(float)*natoms*(nstruct/skip)*3);
  }
  else
  {
  VEL= (float *) malloc (sizeof(float)*lnElt1*(nstruct/skip)*3);  
  }
  ifstream trajvel;
  trajvel.open("inputvel.xyz",ios::in);
  
   for(int i=0;i<(start_frame*(natoms+9));i++)
   {
    getline(trajvel,Dummy);
   }
  
  for(int i=0; i<(nstruct/skip);i=i+1)
  {
    for(int a=0; a<9; a++) //UGLY- careful special for the file from lammps with 10 empty lines - needs fix
    {
      getline(trajvel,Dummy);
    }
    if(Elt1=="all")
    {
    for(int j=0; j< natoms; j++)
    {
      trajvel >> felt >> Dum1 >>Dum2>>Dum3;
        VEL[j*3+i*natoms*3] = Dum1;
        VEL[j*3+i*natoms*3+1] = Dum2;
        VEL[j*3+i*natoms*3+2]= Dum3;
    }
    }
    else if(Elt1=="water")
    {
      int counta=0;
    for(int j=0; j< natoms; j++)
    {
      trajvel >> felt >> Dum1 >>Dum2>>Dum3;
      if(felt=="O"||felt=="H")
      {
        VEL[counta*3+i*lnElt1*3] = Dum1;
        VEL[counta*3+i*lnElt1*3+1] = Dum2;
        VEL[counta*3+i*lnElt1*3+2]= Dum3;
	counta++;
      }
    }      
    }
    else
    {
      int counta=0;
    for(int j=0; j< natoms; j++)
    {
      trajvel >> felt >> Dum1 >>Dum2>>Dum3;
      if(felt==Elt1)
      {
        VEL[counta*3+i*lnElt1*3] = Dum1;
        VEL[counta*3+i*lnElt1*3+1] = Dum2;
        VEL[counta*3+i*lnElt1*3+2]= Dum3;
	counta++;
      }
    }      
    }
    getline(trajvel,Dummy);
    for (int m = 0 ; m < skip-1; m++)
    {
     for(int a=0; a<9; a++) //UGLY- careful special for the file from lammps with 10 empty lines - needs fix
     {
       getline(trajvel,Dummy);
     }      
     for(int j=0;j<natoms;j++)
     {
       getline(trajvel,Dummy);
     }
    }
  }
  trajvel.close();
  }

 
//########################## UNWRAPPER
/*
if( unwrap == "yes")
{
  
  cout << "Unwrapping the trajectory begins"<<endl;
  float lr[3];

for(int i =0; i<(nstruct/skip); i++)
{
 for(int j =0 ; j<lnElt1; j++)
 {
  for(int k =0; k<3;k++)
  {
   Aint[j*3+i*lnElt1*3+k]=0.0;  
  }
  if(i==0)
  {
   for(int k =0; k<3;k++)
   {
    Aint[j*3+i*lnElt1*3+k]=A[j*3+i*lnElt1*3+k];
   }
  }
 }
}

if(cell_type == "orthorhombic")
 {
  for(int i=1; i<(nstruct/skip); i++)
  {//cout<<lnElt1<<endl;
  //cout<<"hi"<<endl;
      for(int j=0;j<lnElt1;j++)
      {//cout<<Elt1<<" ";
          for(int k=0 ;k<3;k++)
          {
          lr[k] = A[j*3+i*lnElt1*3+k] - A[j*3+(i-1)*lnElt1*3+k];
          if(abs(lr[k]) > lattice[k+i*6]/2.0)
          {
               //cout<<"I am working"<<endl;
              if(lr[k] > 0)
              {
              lr[k] = abs(lr[k]) - lattice[k+i*6];
              Aint[j*3+i*lnElt1*3+k] = Aint[j*3+(i-1)*lnElt1*3+k] - abs(lr[k]); 
              }
              else
              {
              lr[k] = abs(lr[k]) - lattice[k+i*6];
              Aint[j*3+i*lnElt1*3+k] = Aint[j*3+(i-1)*lnElt1*3+k] + abs(lr[k]); 
              }     
          }
           else
           {
              Aint[j*3+i*lnElt1*3+k] = Aint[j*3+(i-1)*lnElt1*3+k] + lr[k]; 
           } //cout<<Aint[j*3+i*lnElt1*3+k]<<" ";
          }//cout<<endl;
      }
  }
 }
 else if(cell_type == "monoclinic")
 {
  for(int i=1; i<(nstruct/skip); i++)
  {
      for(int j=0;j<lnElt1;j++)
      {        
        lr[0]=A[j*3+i*lnElt1*3]-A[j*3+(i-1)*lnElt1*3];
        lr[1]=A[j*3+i*lnElt1*3+1]-A[j*3+(i-1)*lnElt1*3+1];
        lr[2]=A[j*3+i*lnElt1*3+2]-A[j*3+(i-1)*lnElt1*3+2];
	
	if(abs(lr[1]) > lattice[1+i*6]/2.0)
        {
        if(lr[1] > 0)
        {
            lr[1] = abs(lr[1]) - lattice[1+i*6];
            Aint[j*3+i*lnElt1*3+1] = Aint[j*3+(i-1)*lnElt1*3+1] - abs(lr[1]);
            lr[0] = abs(lr[0]) - lattice[3+i*6];
            Aint[j*3+i*lnElt1*3+0] = Aint[j*3+(i-1)*lnElt1*3+0] - abs(lr[0]);
	   // Aint[j*3+i*lnElt1*3+0] = Aint[j*3+i*lnElt1*3+0] - lattice[3+i*6];
        }
        else
        {
            lr[1] = abs(lr[1]) - lattice[1+i*6];
            Aint[j*3+i*lnElt1*3+1] = Aint[j*3+(i-1)*lnElt1*3+1] + abs(lr[1]);
            lr[0] = abs(lr[0]) - lattice[3+i*6];
            Aint[j*3+i*lnElt1*3+0] = Aint[j*3+(i-1)*lnElt1*3+0] + abs(lr[0]);
	    //Aint[j*3+i*lnElt1*3+0] = Aint[j*3+i*lnElt1*3+0] + lattice[3+i*6];
        }
        }
        else
        {
        Aint[j*3+i*lnElt1*3+1] = Aint[j*3+(i-1)*lnElt1*3+1] + lr[1];
        }
	
	
        for(int k=0;k<3;k=k+2)
        {
        if(abs(lr[k]) > lattice[k+i*6]/2.0)
        {
            if(lr[k] > 0)
            {
                lr[k]= abs(lr[k]) - lattice[k+i*6];
                Aint[j*3+i*lnElt1*3+k] = Aint[j*3+(i-1)*lnElt1*3+k] - abs(lr[k]);
            }
            else
            {
                lr[k]= abs(lr[k]) - lattice[k+i*6];
                Aint[j*3+i*lnElt1*3+k] = Aint[j*3+(i-1)*lnElt1*3+k] + abs(lr[k]);
            }
        }
        else
        {
            Aint[j*3+i*lnElt1*3+k] = Aint[j*3+(i-1)*lnElt1*3+k] + lr[k];
        }
        }
      }
  }
 }


for(int i =0; i<(nstruct/skip); i++)
{
for(int j =0 ; j<lnElt1; j++)
{
for(int k =0; k<3;k++)
{
A[j*3+i*lnElt1*3+k]=Aint[j*3+i*lnElt1*3+k];
}
}
}

if(unwrapout=="yes")
{

ofstream unwrapped;
unwrapped.open("input_unwrapped.xyz");

for(int i =0; i<(nstruct/skip); i++)
{
unwrapped<<lnElt1<<endl;
unwrapped<<endl;
for(int j =0 ; j<lnElt1; j++)
{unwrapped<<Elt1 <<" ";
for(int k =0; k<3;k++)
{
unwrapped<<Aint[j*3+i*lnElt1*3+k]<<" ";
}
unwrapped<<endl;
}
}
unwrapped.close();
}

free(Aint);
}


//############################# UNWRAPPER FINISH
  */
//in case of same element RDF

  if(Elt1==Elt2)
  {
    B=A;
  }  
  nElt1=lnElt1;
  nElt2=lnElt2;
}




