//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################



#include "cudatools.cuh"
//float minx,maxx,minz,maxz,miny,maxy;
float xrange,yrange,zrange;
float xint,yint,zint;
float *xtick,*ytick,*ztick,*dev_xtick,*dev_ytick,*dev_ztick;
int *density, *dev_density;
void compute_densitygrid()
{
  
  
  ////////////////////Initial setup////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
    if(plot=="int")        { nElt1 = nA_int;}
    else if(plot=="bulk")  { nElt1 = nA_bulk;}
    
cout << type <<" is chosen"<<endl;
cout<<"densitygrid: " << nElt1 <<" "<< Elt1 << " atoms are there in each structure"<<endl;
    printf("%f %f %f %f %f %f\n",minx,maxx,miny,maxy,minz,maxz);
    
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
    xrange=maxx-minx;
    yrange=maxy-miny;
    zrange=maxz-minz;

    xint = xrange/xsplit;
    yint = yrange/ysplit;
    zint = zrange/zsplit;

    xtick=(float*) malloc (sizeof(float)*xsplit*2);
    ytick=(float*) malloc (sizeof(float)*ysplit*2);
    ztick=(float*) malloc (sizeof(float)*zsplit*2);
    for(int i=0;i<xsplit;i++)
    {
      xtick[i*2]=minx+(i*xint);
      xtick[i*2+1]=minx+((i+1)*xint);
    }
    for(int i=0;i<ysplit;i++)
    {
      ytick[i*2]=miny+(i*yint);
      ytick[i*2+1]=miny+((i+1)*yint);
    }
    for(int i=0;i<zsplit;i++)
    {
      ztick[i*2]=minz+(i*zint);
      ztick[i*2+1]=minz+((i+1)*zint);
    }
    
    ofstream plotgnu;
    ofstream plotgnu1;
    plotgnu.open("plot.gnu");
    plotgnu1.open("plot1.gnu");
    
  if(type=="density") 
  {
    
    
     plotgnu1.precision(4);    
     plotgnu1<<"set view map"<<endl\
     <<"r(x) = (x==0 ? 1 : abs(2*x - 0.5)) "<<endl\
     <<"g(x) = (x==0 ? 1 : sin(pi*x))"<<endl\
     <<"b(x) = (x==0 ? 1 : cos( pi*x/2))"<<endl\
     <<"set palette functions r(gray),g(gray),b(gray)"<<endl\
     //<<"set palette defined (0  0.0 0.0 0.5, 1  0.0 0.0 1.0, \
                     2  0.0 0.5 1.0, \
                     3  0.0 1.0 1.0, \
                     4  0.5 1.0 0.5, \
                     5  1.0 1.0 0.0, \
                     6  1.0 0.5 0.0, \
                     7  1.0 0.0 0.0, \
                     8  0.5 0.0 0.0 )"<<endl\
                     
     <<"set pm3d interpolate 20,20"<<endl\
     <<"unset key"<<endl;
     
    plotgnu.precision(4);
    plotgnu<<"unset key"<<endl;
  }
  if(type=="density-top") 
  {
    if(order=="YZX") {
      ztick[0]=mintop;
      ztick[1]=maxtop;      
    }
     if(order=="XYZ") {
      ytick[0]=mintop;
      ytick[1]=maxtop;      
    }
    if(order=="ZXY") {
      xtick[0]=mintop;
      xtick[1]=maxtop;      
    }
    if(choose_atoms != "yes" &&  choose_atoms_from != "yes")
    {
      maxbondist=maxtop;
      minbondist=mintop;
    }
    plotgnu.precision(4);    
    
    plotgnu<<"set view map"<<endl\

     
     <<"r(x) = ( abs(2*x - 0.5)) "<<endl\
     <<"g(x) = ( sin(pi*x))"<<endl\
     <<"b(x) = ( cos( pi*x/2))"<<endl\
     <<"set palette functions r(gray),g(gray),b(gray)"<<endl\

<<"set pm3d interpolate 20,20"<<endl\
//<<"set size ratio @@ratio@@"<<endl\

<<"unset key"<<endl\
<<"#set contour base"<<endl\
<<"#set cntrparam bspline"<<endl\
<<"#set cntrparam levels 5"<<endl;
  }

  FILE *xyztick = fopen("xyztick.dat","wt");
  fprintf(xyztick,"Xsplit");
    for(int i=0;i<xsplit;i++)
    {
      fprintf(xyztick,"%d %f %f \n",i,xtick[i*2],xtick[i*2+1]);
    }
      fprintf(xyztick," \n");
      fprintf(xyztick,"Ysplit");
    for(int i=0;i<ysplit;i++)
    {
      fprintf(xyztick,"%d %f %f \n",i,ytick[i*2],ytick[i*2+1]);    
    }
      fprintf(xyztick," \n");
      fprintf(xyztick,"Zsplit");
    for(int i=0;i<zsplit;i++)
    {
      fprintf(xyztick,"%d %f %f \n",i,ztick[i*2],ztick[i*2+1]);
    }
    fprintf(xyztick," \n");
    int x1= xsplit ;
    int x2= (xsplit+(xsplit*2))/2;
    int x3= (xsplit*2);
    int y1=ysplit;
    int y2=(ysplit*2);
    if(order=="YZX")  {fprintf(xyztick,"%f %f %f %f %f\n",xtick[0],xtick[(xsplit-1)*2+1],ytick[0],ytick[(ysplit-1)*2+1],xtick[(xsplit-1)*2+1]/ytick[(ysplit-1)*2+1]); 
                       plotgnu<<"set size ratio "<<xtick[(xsplit-1)*2+1]/ytick[(ysplit-1)*2+1]<<endl;
		       plotgnu1<<"set size ratio "<<xtick[(xsplit-1)*2+1]/ytick[(ysplit-1)*2+1]<<endl;
                       if(type=="density-top")
                       {
                plotgnu<<"set xrange ["<< xsplit-(xsplit/10)<<":"<< xsplit*2+(xsplit/10)<<"]"<<endl \
                       <<"set yrange ["<< ysplit-(ysplit/10)<<":"<< ysplit*2+(ysplit/10)<<"]"<<endl \
                       <<"set encoding iso_8859_1"<<endl \
		       <<"set xlabel \"Y-direction [\305]\" font 'Verdana,22' offset 0,-1" <<endl \
                       <<"set ylabel \"X-direction [\305]\" font 'Verdana,22' offset -2,0" <<endl \
                       <<"set xtic  font 'Verdana,22' " <<endl \
                       <<"set ytic  font 'Verdana,22' " <<endl \
                        <<"set cbtic  font 'Verdana,22' " <<endl \
                       <<"set cblabel \"Atom density\" font 'Verdana,22' offset 2,0" <<endl \
                       <<"set format x \"%.0f\";set ytics (\"0\" "<< x1<<", \""<<xtick[(xsplit-1)*2+1]/2<<"\" "<<x2<<", \""<< xtick[(xsplit-1)*2+1]<<"\" "<< x3<<" )"<<endl\
                       <<"set format y \"%.0f\";set xtics (\"0\" "<< y1<<", \""<< ytick[(ysplit-1)*2+1]<<"\" "<< y2<<" )"<<endl;}
                       
                       if(type=="density")
                       {
		plotgnu<<"set xrange [ 0 :"<< xtick[(xsplit*2)-1] <<"]"<<endl;
                plotgnu<<"set yrange [ 0 :"<< ytick[(ysplit*2)-1] <<"]"<<endl; 
                plotgnu<<"set xlabel \"Z-direction\" font 'Verdana,22' offset 0,-1" <<endl\
                       <<"set ylabel \"X-direction\" font 'Verdana,22' offset -2,0" <<endl\
                       <<"set xtic  font 'Verdana,20' " <<endl\
                       <<"set ytic  font 'Verdana,20' " <<endl\
                       <<"set xtics 5 "<<endl\
                       <<"set ytics 5"<<endl;
	       plotgnu1<<"set xrange ["<< xsplit-(xsplit/10)<<":"<< xsplit*2+(xsplit/10)<<"]"<<endl\
		       <<"set yrange ["<< ysplit<<":"<< ysplit*2+(ysplit/10)<<"]"<<endl\
                       <<"set xlabel \"Z-direction\" font 'Verdana,22' offset 0,-1" <<endl\
                       <<"set ylabel \"X-direction\" font 'Verdana,22' offset -2,0" <<endl\
                       <<"set xtic  font 'Verdana,22' " <<endl\
                       <<"set ytic  font 'Verdana,22' " <<endl\
                       <<"set cbtic  font 'Verdana,22' " <<endl\
                       <<"set cblabel \"Atom density\" font 'Verdana,22' offset 2,0" <<endl\
                       <<"set format x \"%.0f\";set ytics (\"0\" "<< x1<<", \""<<xtick[(xsplit-1)*2+1]/2<<"\" "<<x2<<", \""<< xtick[(xsplit-1)*2+1]<<"\" "<< x3<<" )"<<endl\
                       <<"set format y \"%.0f\";set xtics (\"0\" "<< y1<<",  \""<< ytick[(ysplit-1)*2+1]<<"\" "<< y2<<" )"<<endl\
		       <<"unset cblabel; unset colorbox"<<endl;
		      }
                      }
                       
    else if(order=="XYZ") { fprintf(xyztick,"%f %f \n",ztick[0],ztick[(zsplit-1)*2+1],xtick[0],xtick[(xsplit-1)*2+1],ztick[(zsplit-1)*2+1]/xtick[(xsplit-1)*2+1]);
                       plotgnu<<"set size ratio "<<ztick[(zsplit-1)*2+1]/xtick[(xsplit-1)*2+1]<<endl;
		       plotgnu1<<"set size ratio "<<ztick[(zsplit-1)*2+1]/xtick[(xsplit-1)*2+1]<<endl;
                       if(type=="density-top")
                       {
                plotgnu<<"set xrange ["<< xsplit-(xsplit/10)<<":"<< xsplit*2+(xsplit/10)<<"]"<<endl\
                       <<"set yrange ["<< zsplit-(zsplit/10)<<":"<< zsplit*2+(zsplit/10)<<"]"<<endl\
                       <<"set encoding iso_8859_1"<<endl \
		       <<"set xlabel \"Y-direction [\305]\" font 'Verdana,22' offset 0,-1" <<endl \
                       <<"set ylabel \"X-direction [\305]\" font 'Verdana,22' offset -2,0" <<endl \
                       <<"set xtic  font 'Verdana,22' " <<endl\
                       <<"set ytic  font 'Verdana,22' " <<endl\
                        <<"set cbtic  font 'Verdana,22' " <<endl\
                       <<"set cblabel \"Atom density\" font 'Verdana,22' offset 2,0" <<endl\
                       <<"set format x \"%.0f\";set ytics (\"0\" "<< zsplit<<", \""<<ztick[(zsplit-1)*2+1]/2<<"\" "<<(zsplit+(zsplit*2))/2<<", \""<< ztick[(zsplit-1)*2+1]<<"\" "<< (zsplit*2)<<" )"<<endl\
                       <<"set format y \"%.0f\";set xtics (\"0\" "<< x1<<",  \""<< xtick[(xsplit-1)*2+1]<<"\" "<< (xsplit*2)<<" )"<<endl;
			 
		      }
                       
                       
                       if(type=="density")
                       {
		plotgnu<<"set xrange [0:"<< xtick[(xsplit*2)-1] <<"]"<<endl\
                       <<"set yrange [0:"<< ztick[(zsplit*2)-1] <<"]"<<endl\
                       <<"set xlabel \"X-direction\" font 'Verdana,22' offset 0,-1" <<endl\
                       <<"set ylabel \"Z-direction\" font 'Verdana,22' offset -2,0" <<endl\
                       <<"set xtic  font 'Verdana,20' " <<endl\
                       <<"set ytic  font 'Verdana,20' " <<endl\
                       <<"set xtics 5 "<<endl\
                       <<"set ytics 5"<<endl;
		plotgnu1<<"set xrange ["<< xsplit-(xsplit/10)<<":"<< xsplit*2+(xsplit/10)<<"]"<<endl\
		<<"set yrange ["<< zsplit<<":"<< zsplit*2+(zsplit/10)<<"]"<<endl\
                       <<"set xlabel \"X-direction\" font 'Verdana,22' offset 0,-1" <<endl\
                       <<"set ylabel \"Z-direction\" font 'Verdana,22' offset -2,0" <<endl\
                       <<"set xtic  font 'Verdana,20'" <<endl\
                       <<"set ytic  font 'Verdana,20' " <<endl\
                        <<"set cbtic  font 'Verdana,20' " <<endl\
                       <<"set cblabel \"Atom density\" font 'Verdana,20' offset 2,0" <<endl\
                       <<"set format x \"%.0f\";set ytics (\"0\" "<< zsplit<<", \""<<ztick[(zsplit-1)*2+1]/2<<"\" "<<(zsplit+(zsplit*2))/2<<", \""<< ztick[(zsplit-1)*2+1]<<"\" "<< (zsplit*2)<<" )"<<endl\
                       <<"set format y \"%.0f\";set xtics (\"0\" "<< x1<<",  \""<< xtick[(xsplit-1)*2+1]<<"\" "<< (xsplit*2)<<" )"<<endl\
		       <<"unset cblabel; unset colorbox"<<endl; 
		      }
                       }
    
    else if(order=="ZXY") { fprintf(xyztick,"%f %f \n",ytick[0],ytick[(ysplit-1)*2+1],ztick[0],ztick[(zsplit-1)*2+1],ytick[(ysplit-1)*2+1]/ztick[(zsplit-1)*2+1]);
                       plotgnu<<"set size ratio "<<ytick[(ysplit-1)*2+1]/ztick[(zsplit-1)*2+1]<<endl;
		       plotgnu1<<"set size ratio "<<ytick[(ysplit-1)*2+1]/ztick[(zsplit-1)*2+1]<<endl;
                       if(type=="density-top")
                       {
                       plotgnu<<"set xrange ["<< zsplit-(zsplit/10)<<":"<< zsplit*2+(zsplit/10)<<"]"<<endl\
                       <<"set yrange ["<< ysplit-(ysplit/10)<<":"<< ysplit*2+(ysplit/10)<<"]"<<endl\
                       <<"set encoding iso_8859_1"<<endl \
		       <<"set xlabel \"Y-direction [\305]\" font 'Verdana,22' offset 0,-1" <<endl \
                       <<"set ylabel \"X-direction [\305]\" font 'Verdana,22' offset -2,0" <<endl \
                       <<"set xtic  font 'Verdana,22' " <<endl\
                       <<"set ytic  font 'Verdana,22' " <<endl\
                       <<"set cbtic  font 'Verdana,22' " <<endl\
                       <<"set cblabel \"Atom density\" font 'Verdana,22' offset 2,0" <<endl\
                       <<"set format x \"%.0f\";set ytics (\"0\" "<< ysplit<<", \""<<ytick[(ysplit-1)*2+1]/2<<"\" "<<(ysplit+(ysplit*2))/2<<", \""<< ytick[(ysplit-1)*2+1]<<"\" "<< (ysplit*2)<<" )"<<endl\
                       <<"set format y \"%.0f\";set xtics (\"0\" "<< zsplit<<",  \""<< ztick[(zsplit-1)*2+1]<<"\" "<< (zsplit*2)<<" )"<<endl;}
                       
                       if(type=="density")
                       {
	        plotgnu<<"set xrange [0:"<< ztick[(zsplit*2)-1] <<"]"<<endl\
                       <<"set yrange [0:"<< ytick[(ysplit*2)-1] <<"]"<<endl\
                       <<"set xlabel \"Y-direction\" font 'Verdana,22' offset 0,-1" <<endl\
                       <<"set ylabel \"Z-direction\" font 'Verdana,22' offset -2,0" <<endl\
                       <<"set xtic  font 'Verdana,20' " <<endl\
                       <<"set ytic  font 'Verdana,20' " <<endl\
                       <<"set xtics 5 "<<endl\
                       <<"set ytics 5"<<endl;
		plotgnu1<<"set xrange ["<< zsplit-(zsplit/10)<<":"<< zsplit*2+(zsplit/10)<<"]"<<endl\
		<<"set yrange ["<< ysplit<<":"<< ysplit*2+(ysplit/10)<<"]"<<endl\
                       <<"set xlabel \"Y-direction\" font 'Verdana,22' offset 0,-1" <<endl\
                       <<"set ylabel \"Z-direction\" font 'Verdana,22' offset -2,0" <<endl\
                       <<"set xtic  font 'Verdana,22' " <<endl\
                       <<"set ytic  font 'Verdana,22' " <<endl\
                       <<"set cbtic  font 'Verdana,22' " <<endl\
                       <<"set cblabel \"Atom density\" font 'Verdana,22' offset 2,0" <<endl\
                       <<"set format x \"%.0f\";set ytics (\"0\" "<< ysplit<<", \""<<ytick[(ysplit-1)*2+1]/2<<"\" "<<(ysplit+(ysplit*2))/2<<", \""<< ytick[(ysplit-1)*2+1]<<"\" "<< (ysplit*2)<<" )"<<endl\
                       <<"set format y \"%.0f\";set xtics (\"0\" "<< zsplit<<",  \""<< ztick[(zsplit-1)*2+1]<<"\" "<< (zsplit*2)<<" )"<<endl\
		       <<"unset cblabel; unset colorbox"<<endl;
		      }
                       }
fclose(xyztick);
plotgnu<<"set terminal postscript solid eps enhanced color font 'Verdana,18'"<<endl\
       <<"set output 'Plot.eps'"<<endl;

       if(type=="density-top")
       {
       plotgnu<<"sp 'extended.out' i 0 matrix w pm3d , 'teter.dat' u ($1)+49.5:($2)+49.5:(1) w p pt 65 ps 2 lc rgb \"black\" lw 3"<<endl;
       }
       else if(type=="density")
	 {
	  
	  plotgnu1<<"set terminal postscript solid eps enhanced color font 'Verdana,18'"<<endl\
          <<"set output 'Plot1.eps'"<<endl;
	  plotgnu1<<"sp 'extended.out' i 0 matrix w pm3d"<<endl;
	 }
plotgnu.close();
plotgnu1.close();
int *densityf, *dev_densityf;
    density=(int *) malloc (sizeof(int)*xsplit*ysplit*zsplit);
    densityf=(int *) malloc (sizeof(int)*xsplit*ysplit*zsplit);
    for(int i=0;i<xsplit;i++)
    {
      for(int j=0;j<ysplit;j++)
      {
        for(int k=0;k<zsplit;k++)
        {
          density[k+j*zsplit+i*ysplit*zsplit]=0;
	  densityf[k+j*zsplit+i*ysplit*zsplit]=0;
        }
      }
    }
///*

///////////////////////////////////////////Initial setup done/////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////Special atom selection between 2 distances/////////////////////////////////////////////////////////////////////////
float surf,avgsurf=0.0;
int dirn;
  if(choose_dirn=="z") dirn=2;
  else if(choose_dirn=="x") dirn=0;
  else if(choose_dirn=="y") dirn=1;
  
if(choose_atoms=="yes")
{
cout << Elt1 <<" within "<< minbondist << " and " << maxbondist<<" angstrom of "<< metal_species<<endl;
  
 
  
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
    cudaMalloc((void **)&dev_B,sizeof(float)*count_metal*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_exch,sizeof(int)*nElt1*(nstruct/skip));
    cudaMalloc((void **)&dev_lattice,sizeof(float)*6*(nstruct/skip));
//Copy data from host to device

    cudaMemcpy(dev_A,A,sizeof(float)*nElt1*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_B,METAL,sizeof(float)*count_metal*(nstruct/skip)*3,cudaMemcpyHostToDevice);

//dim3 struct to define elements of the execution configuration


    dim3 dimBlocke(10,10,10);
    dim3 dimGride((nElt1+dimBlocke.x-1)/dimBlocke.x,(count_metal+dimBlocke.y-1)/dimBlocke.y,((nstruct/skip)+dimBlocke.z-1)/dimBlocke.z);

//Cuda kernal execution for distance matrix with CUDA timing API commands

    //For surface metal atoms
    cudaMemcpy(dev_exch,exch,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice);
    cudaMemcpy(dev_lattice,lattice,sizeof(int)*6*(nstruct/skip),cudaMemcpyHostToDevice);
    if (cell_type == "orthorhombic") 
    {
      if(strict_dirn == "yes") {covermat<<<dimGride,dimBlocke>>>(dev_A,dev_exch,minbondist,maxbondist,nElt1,(nstruct/skip),dirn,avgsurf,mintop,maxtop);}
      else {covermat<<<dimGride,dimBlocke>>>(dev_A,dev_B,dev_exch,minbondist,maxbondist,nElt1,count_metal,(nstruct/skip),dev_lattice,dirn,avgsurf,mintop,maxtop);}
    }
    else if (cell_type == "monoclinic")
    {
      if(strict_dirn == "yes") {covermat<<<dimGride,dimBlocke>>>(dev_A,dev_exch,minbondist,maxbondist,nElt1,(nstruct/skip),dirn,avgsurf,mintop,maxtop);} //no lattice needed thus use covermat - no problem
      else {covermatmono<<<dimGride,dimBlocke>>>(dev_A,dev_B,dev_exch,minbondist,maxbondist,nElt1,count_metal,(nstruct/skip),dev_lattice,dirn,avgsurf,mintop,maxtop);}
    }
    
    
    cudaMemcpy(exch,dev_exch,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyDeviceToHost);
    
    
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
          if (exch[i+l*nElt1]>0) {exch[i+l*nElt1]=1; Aintnum[l]++;}
          else if (exch[i+l*nElt1]==0) Abulknum[l]++;
	  //exchout<<i<<" "<<exch[i+l*nElt1]<<endl;
        }
        //exchout<<endl;
      }
      
      //exchout.close();
     cudaMemcpy(dev_exch,exch,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice); 
     
     ofstream coverage;
     coverage.open("coverage.data");
     for(int l=0; l<nstruct/skip; l++)
      {
	coverage <<l<<" "<< float(Aintnum[l]) <<" "<<float(Aintnum[l])/float(surfatom)<<endl;
      }
     coverage.close();
     
}          
   
///////////////////////////////////////////////special selection done - stored in exch///////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////Special atom selection between 2 distances from a particlar index/////////////////////////////////////////////////////////////////////////

if(choose_atoms_from=="yes")
{
cout << Elt1 <<" within "<< minbondist << " and " << maxbondist<<" angstrom of "<< check_index<<endl;
  
  float *index_xyz;
  index_xyz = (float *) malloc (sizeof(float)*(nstruct/skip)*3);
  for(int i=0 ; i<(nstruct/skip); i++)
  {
  for(int j=0 ; j<3; j++)
  {
    index_xyz[i*3+j]=0.0;
  }
  }
  
  //careful!! directly we choose metal array!! not good!! but works for now!!
  
  for(int i=0 ; i<(nstruct/skip); i++)
  {
  for(int j=0 ; j<3; j++)
  {
    index_xyz[i*3+j]=METAL[check_index*3+i*count_metal*3+j];
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
    cudaMalloc((void **)&dev_exch,sizeof(int)*nElt1*(nstruct/skip));
    cudaMalloc((void **)&dev_lattice,sizeof(float)*6*(nstruct/skip));
//Copy data from host to device

    cudaMemcpy(dev_A,A,sizeof(float)*nElt1*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_B,index_xyz,sizeof(float)*(nstruct/skip)*3,cudaMemcpyHostToDevice);

//dim3 struct to define elements of the execution configuration

    
    dim3 dimBlocke(32,1,32);
    dim3 dimGride((nElt1+dimBlocke.x-1)/dimBlocke.x,1,((nstruct/skip)+dimBlocke.z-1)/dimBlocke.z);

//Cuda kernal execution for distance matrix with CUDA timing API commands

    //For surface metal atoms
    cudaMemcpy(dev_exch,exch,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice);
    cudaMemcpy(dev_lattice,lattice,sizeof(int)*6*(nstruct/skip),cudaMemcpyHostToDevice);
    if (cell_type == "orthorhombic") 
    {
      covermat<<<dimGride,dimBlocke>>>(dev_A,dev_B,dev_exch,minbondist,maxbondist,check_index,dev_lattice,nElt1,(nstruct/skip));
    }
    else if (cell_type == "monoclinic")
    {
      covermatmono<<<dimGride,dimBlocke>>>(dev_A,dev_B,dev_exch,minbondist,maxbondist,check_index,dev_lattice,nElt1,(nstruct/skip));
    }
    
    
    cudaMemcpy(exch,dev_exch,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyDeviceToHost);
    
    
    cudaFree(dev_A);
    cudaFree(dev_B);
    cudaFree(dev_lattice);
    
    ofstream exchout;
    exchout.open("exchout.data");
    
    
    for(int l=0; l<(nstruct/skip); l++)
      {
	//exchout<<l<<" "<<index_xyz[l*3]<<" "<<index_xyz[l*3+1]<<" "<<index_xyz[l*3+2]<<endl;
      for(int i=0; i< nElt1; i++)
        {
	  //exchout<<i<<" "<<exch[i+l*nElt1]<<endl;
          if (exch[i+l*nElt1]>0) {exch[i+l*nElt1]=1; Aintnum[l]++;}
          else if (exch[i+l*nElt1]==0) Abulknum[l]++;
	  //exchout<<i<<" "<<exch[i+l*nElt1]<<endl;
        }
        exchout<<l<<" "<<Aintnum[l]<<" "<<Abulknum[l]<<endl;
      }
      
      exchout.close();
      
     cudaMemcpy(dev_exch,exch,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice); 

     
}          
   
///////////////////////////////////////////////special selection done - stored in exch///////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////Density calculation starts/////////////////////////////////////////////////////////////////////////////////////////////

if(density_plot=="yes")
{

    printf("Start of cuda calculation\n");
    cudaMalloc((void **)&dev_A,sizeof(float)*(nstruct/skip)*nElt1*3);
    cudaMalloc((void **)&dev_density,sizeof(int)*xsplit*ysplit*zsplit);
    cudaMalloc((void **)&dev_xtick, sizeof(float)*xsplit*2);
    cudaMalloc((void **)&dev_ytick, sizeof(float)*ysplit*2);
    cudaMalloc((void **)&dev_ztick, sizeof(float)*zsplit*2);
    if(plot=="int")
    {
      cudaMemcpy(dev_A,A_int,sizeof(float)*(nstruct/skip)*nA_int*3,cudaMemcpyHostToDevice);
    }
    else if(plot=="bulk")
    {
      cudaMemcpy(dev_A,A_bulk,sizeof(float)*(nstruct/skip)*nA_bulk*3,cudaMemcpyHostToDevice);
    }
    else if(plot=="all")
    {
      cudaMemcpy(dev_A,A,sizeof(float)*(nstruct/skip)*nElt1*3,cudaMemcpyHostToDevice);
    }
    cudaMemcpy(dev_density,density,sizeof(int)*xsplit*ysplit*zsplit,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_xtick,xtick,sizeof(float)*xsplit*2,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_ytick,ytick,sizeof(float)*ysplit*2,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_ztick,ztick,sizeof(float)*zsplit*2,cudaMemcpyHostToDevice);
    
    printf("%ld B needed\n",(((nstruct/skip)*nElt1*3+xsplit*ysplit*zsplit+xsplit*2+ysplit*2+zsplit*2)*4));
    dim3 dimBlock(32,32,1);
    //dim3 dimBlock(1,1,1);
    dim3 dimGrid(((nstruct/skip)+dimBlock.x-1)/dimBlock.x,(nElt1+dimBlock.y-1)/dimBlock.y,1);
    //dim3 dimGrid(nstruct,nElt1,1);
    if(choose_atoms=="yes")
    {
    density_calc<<<dimGrid,dimBlock>>>(dev_A,dev_density,dev_exch,(nstruct/skip),nElt1,xsplit,ysplit,zsplit,dev_xtick,dev_ytick,dev_ztick);
    cudaFree(dev_exch);
    }
    else if(choose_atoms_from=="yes")
    {
    density_calc<<<dimGrid,dimBlock>>>(dev_A,dev_density,dev_exch,(nstruct/skip),nElt1,xsplit,ysplit,zsplit,dev_xtick,dev_ytick,dev_ztick);
    cudaFree(dev_exch);
    }
    else
    {
      cout<<"I am computing"<<endl;
    density_calc<<<dimGrid,dimBlock>>>(dev_A,dev_density,(nstruct/skip),nElt1,xsplit,ysplit,zsplit,dev_xtick,dev_ytick,dev_ztick);      
    }
    //cudaMemcpy(density,dev_density,sizeof(int)*xsplit*ysplit*zsplit,cudaMemcpyDeviceToHost);
     cudaFree(dev_A);
}
/////////////////////////////////////////////////////Density calculation done//////////////////////////////////////////////////////////////////////////////////






////////////////////////////////////////////////////Dipole orientation computation/////////////////////////////////////////////////////////////////////////////     
    
   if(dipole_orient=="yes" && type=="density-top")
   {
    
    int *h2onum;
  int *h2oneigh;
  float *h2oxyz;
  float *h2oneighdist;
  float *h2odistmat, *dev_h2odistmat;
  
  int OHsamples=20;
  
    h2onum = (int *) malloc (sizeof(int)*nElt1*(nstruct/skip));
    h2oneigh = (int *) malloc (sizeof(int)*OHsamples*nElt1*(nstruct/skip));
    h2oneighdist = (float *) malloc (sizeof(float)*OHsamples*nElt1*(nstruct/skip));
    h2oxyz = (float *) malloc (sizeof(float)*15*nElt1*(nstruct/skip));  //Oxyz,H1xyz,H2xyz,Hmidpoint,dipole alignment with surface normal
    h2odistmat = (float *) malloc (sizeof(float)*nElt1*nElt2);
    
  for(int i=0;i<nstruct/skip;i++)
  {
    for(int j=0; j<nElt1; j++)
    {
      h2onum[j+i*nElt1] =0;
      for(int k =0 ;k<OHsamples;k++)
      {
      h2oneigh[j*OHsamples+i*nElt1*OHsamples+k]=0;
      h2oneighdist[j*OHsamples+i*nElt1*OHsamples+k]=0.0;
      }
      for(int k =0 ;k<15;k++)
      {
      h2oxyz[j*15+i*nElt1*15+k]=0.0;
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
    
     
    //ofstream neigh;
   // neigh.open("H2O.dat");
   
    
for(int i = 0; i<(nstruct/skip); i++)
{
  
  //neigh << i<<endl<<endl;
  
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
          Hbondmat<<<dimGrid,dimBlock>>>(dev_A,dev_B,dev_h2odistmat,bondistOH,nElt1,nElt2,dev_lattice,i);
    }
    else if(cell_type == "monoclinic")
    {
          Hbondmatmono<<<dimGrid,dimBlock>>>(dev_A,dev_B,dev_h2odistmat,bondistOH,nElt1,nElt2,dev_lattice,i);      
    }

cudaMemcpy(h2odistmat,dev_h2odistmat,sizeof(float)*nElt1*nElt2,cudaMemcpyDeviceToHost);
   
  
  for(int i1=0;i1<nElt1;i1++) //O
  {
    int k=0;
    //neigh << i1<<" ";
    for(int j1=0; j1<nElt2; j1++) //H
    {
      //neigh <<j1<<" "<<h2odistmat[j1+i1*nElt2]<<" ";
      if(h2odistmat[j1+i1*nElt2] !=0.0)
      {
	h2onum[i1+i*nElt1]++;
	h2oneigh[i1*OHsamples+i*nElt1*OHsamples+k] = j1;
	h2oneighdist[i1*OHsamples+i*nElt1*OHsamples+k] = h2odistmat[j1+i1*nElt2];
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
  
    //sorting needs to be done!!!!

cout<<"Swapping of bonds based on distance "<<endl;

float swapa,swapb;
    for(int i=0;i<(nstruct/skip);i++)
     {
      for(int j=0; j<nElt1; j++)
       {
	 
	for(int swf=0; swf<h2onum[j+i*nElt1]-1;swf++)
	{
	for(int sws=swf+1; sws<h2onum[j+i*nElt1];sws++)
	{
	if(h2oneighdist[j*OHsamples+i*nElt1*OHsamples+swf] > h2oneighdist[j*OHsamples+i*nElt1*OHsamples+sws] && sws !=swf && h2oneighdist[j*OHsamples+i*nElt1*OHsamples+sws] !=0.0 && h2oneighdist[j*OHsamples+i*nElt1*OHsamples+swf] !=0.0) 
	{
	   swapa = h2oneigh[j*OHsamples+i*nElt1*OHsamples+swf]; 
	    h2oneigh[j*OHsamples+i*nElt1*OHsamples+swf]=h2oneigh[j*OHsamples+i*nElt1*OHsamples+sws];
	    h2oneigh[j*OHsamples+i*nElt1*OHsamples+sws]=swapa;
	   swapb = h2oneighdist[j*OHsamples+i*nElt1*OHsamples+swf]; 
	    h2oneighdist[j*OHsamples+i*nElt1*OHsamples+swf]=h2oneighdist[j*OHsamples+i*nElt1*OHsamples+sws];
	    h2oneighdist[j*OHsamples+i*nElt1*OHsamples+sws]=swapb;
	}
	}
	}	
       }
     }	
   
cout<<" done !!!"<<endl;
//H2onum has number of H around O within 1.2 AA
//swap H1 short OH and H2, long

 //cout<<"Write Hxyz.dat"<<endl;

  //ofstream watercount;
  //watercount.open("Hxyz.dat");
  
  
  float *angdist,*angdistH1,*angdistH2;
  
  angdist = (float *) malloc (sizeof(float)*180);
  angdistH1 = (float *) malloc (sizeof(float)*180);
  angdistH2 = (float *) malloc (sizeof(float)*180);
  
  for(int i=0; i<=180;i++)
  {
    angdist[i]=0.0;
    angdistH1[i]=0.0;
    angdistH2[i]=0.0;
  }
  
  float val;
  
  
  ofstream OHzdist;
  OHzdist.open("OHzdist.dat");
  ofstream dipOHdist;
  dipOHdist.open("dipOHdist.dat");
for(int i=0; i<(nstruct/skip); i++)
{
 // cout<<i<<endl;
 for(int j=0; j<nElt1;j++)
 {
    if(order=="YZX") {
      val=A[j*3+i*nElt1*3+2];    
    }
    if(order=="XYZ") {
      val=A[j*3+i*nElt1*3+1];      
    }
    if(order=="ZXY") {
      val=A[j*3+i*nElt1*3];      
    }
   if(h2onum[j+i*nElt1]==2 && exch[j+i*nElt1]==1 &&val<maxtop && val>mintop)
   //if(h2onum[j+i*nElt1]==2)
   {
        //watercount << i <<" "<<j<<" ";
     for(int k=0; k<3;k++)
     {
     h2oxyz[j*15+i*nElt1*15+k]= A[j*3+i*nElt1*3+k];
     //watercount << h2oxyz[j*15+i*nElt1*15+k] <<" ";
     }
     for(int k=0; k<3;k++)
     {
     h2oxyz[j*15+i*nElt1*15+3+k]= B[h2oneigh[j*OHsamples+i*nElt1*OHsamples]*3+i*nElt2*3+k];
     //watercount << h2oxyz[j*9+i*nElt1*9+3+k] <<" ";
     }
     for(int k=0; k<3;k++)
     {
     h2oxyz[j*15+i*nElt1*15+6+k]= B[h2oneigh[j*OHsamples+i*nElt1*OHsamples+1]*3+i*nElt2*3+k];
     //watercount << h2oxyz[j*9+i*nElt1*9+6+k] <<" ";
     }
     OHzdist <<A[j*3+i*nElt1*3+2]-avgsurf<<" "<<B[h2oneigh[j*OHsamples+i*nElt1*OHsamples]*3+i*nElt2*3+2]-avgsurf<<endl;
     OHzdist <<A[j*3+i*nElt1*3+2]-avgsurf<<" "<<B[h2oneigh[j*OHsamples+i*nElt1*OHsamples+1]*3+i*nElt2*3+2]-avgsurf<<endl;
     float chk1,chk2,chk3,chk4,chk5,chk6;
     float a1,b1,c1,d1,e1,f1,g1,h1,i1;
     float latx,laty,latz,latxy,latxz,latyz;
     //check_together(h2oxyz[j*15+i*nElt1*15],h2oxyz[j*15+i*nElt1*15+1],h2oxyz[j*15+i*nElt1*15+2],\
                    h2oxyz[j*15+i*nElt1*15+3],h2oxyz[j*15+i*nElt1*15+4],h2oxyz[j*15+i*nElt1*15+5],\
                    h2oxyz[j*15+i*nElt1*15+6],h2oxyz[j*15+i*nElt1*15+7],h2oxyz[j*15+i*nElt1*15+8],\
                    lattice[i*6],lattice[i*6+1],lattice[i*6+2],lattice[i*6+3],lattice[i*9+4],lattice[i*6+5]);
     
     a1=h2oxyz[j*15+i*nElt1*15];
     b1=h2oxyz[j*15+i*nElt1*15+1];
     c1=h2oxyz[j*15+i*nElt1*15+2];
     d1=h2oxyz[j*15+i*nElt1*15+3];
     e1=h2oxyz[j*15+i*nElt1*15+4];
     f1=h2oxyz[j*15+i*nElt1*15+5];
     g1=h2oxyz[j*15+i*nElt1*15+6];
     h1=h2oxyz[j*15+i*nElt1*15+7];
     i1=h2oxyz[j*15+i*nElt1*15+8];
     latx=lattice[i*6];
     laty=lattice[i*6+1];
     latz=lattice[i*6+2];
     latxy=lattice[i*6+3];
     latxz=lattice[i*6+4];
     latyz=lattice[i*6+5];
     
chk1= (a1-d1);
chk2= (b1-e1);
chk3= (c1-f1);
chk4= (a1-g1);
chk5= (b1-h1);
chk6= (c1-i1);

if(cell_type == "orthorhombic")
{
if(fabs(chk1) > latx/2.0) { if (chk1 >0) d1=d1+latx; else d1=d1-latx;}
if(fabs(chk2) > laty/2.0) { if (chk2 >0) e1=e1+laty; else e1=e1-laty;}
if(fabs(chk3) > latz/2.0) { if (chk3 >0) f1=f1+latz; else f1=f1-latz;}
if(fabs(chk4) > latx/2.0) { if (chk4 >0) g1=g1+latx; else g1=g1-latx;}
if(fabs(chk5) > laty/2.0) { if (chk5 >0) h1=h1+laty; else h1=h1-laty;}
if(fabs(chk6) > latz/2.0) { if (chk6 >0) i1=i1+latz; else i1=i1-latz;}
}

else if(cell_type == "monoclinic")
{
if(fabs(chk2) > laty/2.0) { if (chk2 >0) {e1=e1+laty; d1=d1+latxy;} else {e1=e1-laty; d1=d1-latxy;}}
chk1= (a1-d1);
if(fabs(chk1) > latx/2.0) { if (chk1 >0) d1=d1+latx; else d1=d1-latx;}
if(fabs(chk3) > latz/2.0) { if (chk3 >0) f1=f1+latz; else f1=f1-latz;}
if(fabs(chk5) > laty/2.0) { if (chk5 >0) {h1=h1+laty; g1=g1+latxy;} else {h1=h1-laty; g1=g1-latxy;}}
chk4= (a1-g1);
if(fabs(chk4) > latx/2.0) { if (chk4 >0) g1=g1+latx; else g1=g1-latx;}
if(fabs(chk6) > latz/2.0) { if (chk6 >0) i1=i1+latz; else i1=i1-latz;}
}
     
     h2oxyz[j*15+i*nElt1*15+3]=d1;
     h2oxyz[j*15+i*nElt1*15+4]=e1;
     h2oxyz[j*15+i*nElt1*15+5]=f1;
     h2oxyz[j*15+i*nElt1*15+6]=g1;
     h2oxyz[j*15+i*nElt1*15+7]=h1;
     h2oxyz[j*15+i*nElt1*15+8]=i1;
     
     for(int k=0; k<3;k++)
     {
     h2oxyz[j*15+i*nElt1*15+9+k]=(h2oxyz[j*15+i*nElt1*15+3+k]+h2oxyz[j*15+i*nElt1*15+6+k])/2.0;
     //watercount << h2oxyz[j*15+i*nElt1*15+9+k] <<" ";
     }
     //////////////////////////////angle between dipole and z axis 0,0,1(specific - be careful if using other directions//////////////////////////////////////////////////
     h2oxyz[j*15+i*nElt1*15+12]=angle(0,0,0,0,0,1,h2oxyz[j*15+i*nElt1*15],h2oxyz[j*15+i*nElt1*15+1],h2oxyz[j*15+i*nElt1*15+2],h2oxyz[j*15+i*nElt1*15+9],h2oxyz[j*15+i*nElt1*15+9+1],h2oxyz[j*15+i*nElt1*15+9+2]);
     h2oxyz[j*15+i*nElt1*15+13]=angle(0,0,0,0,0,1,h2oxyz[j*15+i*nElt1*15],h2oxyz[j*15+i*nElt1*15+1],h2oxyz[j*15+i*nElt1*15+2],h2oxyz[j*15+i*nElt1*15+3],h2oxyz[j*15+i*nElt1*15+4],h2oxyz[j*15+i*nElt1*15+5]);
     h2oxyz[j*15+i*nElt1*15+14]=angle(0,0,0,0,0,1,h2oxyz[j*15+i*nElt1*15],h2oxyz[j*15+i*nElt1*15+1],h2oxyz[j*15+i*nElt1*15+2],h2oxyz[j*15+i*nElt1*15+6],h2oxyz[j*15+i*nElt1*15+7],h2oxyz[j*15+i*nElt1*15+8]);
     //watercount << h2oxyz[j*15+i*nElt1*15+12]*(180.0/3.141);
     angdist[int(ceil(h2oxyz[j*15+i*nElt1*15+12]*(180.0/3.14159)))]++;
     angdistH1[int(ceil(h2oxyz[j*15+i*nElt1*15+13]*(180.0/3.141)))]++;
     angdistH2[int(ceil(h2oxyz[j*15+i*nElt1*15+14]*(180.0/3.141)))]++;
     //watercount<<endl;
     dipOHdist << j<<" "<<h2oxyz[j*15+i*nElt1*15+12]*(180.0/3.14159)<<" "<<h2oxyz[j*15+i*nElt1*15+13]*(180.0/3.141)<<" "<<\
     h2oxyz[j*15+i*nElt1*15+14]*(180.0/3.141)<<endl;
   }
   
  }
}
  //watercount.close();
  OHzdist.close();
  dipOHdist.close();
  float totalangdist=0.0;
  for(int i=0; i<=180;i++)
  {
    totalangdist +=angdist[i];
  }
  ofstream angldist;
  angldist.open("Dipole-orient.dat");
  angldist <<"T num num1 num2"<<endl;
   for(int i=0; i<=180;i++)
  {
    //angldist << i<<" "<<angdist[i]/(totalangdist)<<endl;
    angldist << i<<" "<<angdist[i]/float(nstruct/skip)<<" "<< angdistH1[i]/float(nstruct/skip)<<" "<< angdistH2[i]/float(nstruct/skip)<<endl;
 
  }
  angldist.close();
  free(angdist);
  free(angdistH1);
  free(angdistH2);
  free(h2onum);
  free(h2oneigh);
  free(h2oneighdist);
  free(h2oxyz);
  free(h2odistmat);
   }
    //distribution of hydrogen atoms per chosen range
    int nlines;
    float DOG;
    float *OHzdistmat;
    int *OHzhist;
    system("wc -l < OHzdist.dat > nlines");
    ifstream nline;
    nline.open("nlines");
    nline>>nlines;
    nline.close();
    system("rm nlines");
    OHzdistmat = (float*) malloc (sizeof(float)*nlines);
    ifstream hist;
    hist.open("OHzdist.dat");
    for(int i=0;i<nlines;i++)
    {
            hist >> DOG >>OHzdistmat[i];
    }
    hist.close();
    //cout<<"hello"<<endl;
    float Hmax,Hint;
    float *Htick;
    Hmax=5.0;
    //Hmax = maximum(OHzdistmat,nlines,0);
     //cout<<"okay"<<endl;
    OHzhist=(int *) malloc (sizeof(int)*100);
    Htick=(float*) malloc (sizeof(float)*100*2);
    //cout<<Hmax<<endl;
    Hint = Hmax/100.0;
    
    for(int i=0;i<100;i++)
    {
      OHzhist[i]=0;
      Htick[i*2]=(i*Hint);
      Htick[i*2+1]=((i+1)*Hint);
    }
    
   
    for(int i=0;i<100;i++)
    {
      for(int j=0;j<nlines;j++)
    {
            if(OHzdistmat[j]>=Htick[i*2] && OHzdistmat[j]<Htick[i*2+1])
	    {
	      OHzhist[i]++;
	    }
    }
    }
    float cellvolume;
    
    if(cell_type=="orthorhombic")
    {
      cellvolume = xvec*yvec*Hint;
    }
    else if(cell_type=="monoclinic")
    {
      cellvolume= xvec*xvec*zint*0.866025403;
    }
    ofstream angldist;
  angldist.open("OHzhist.dat");
  angldist <<"Z num"<<endl;
   for(int i=0; i<100;i++)
  {
    angldist << Htick[i*2+1]<<" "<<float(OHzhist[i])/float((nstruct/skip)*cellvolume)<<endl;
  }
  angldist.close();

  free(OHzdistmat);
  free(OHzhist);
  free(Htick);
////////////////////////////////////////////////////end of dipole orientation calculation /////////////////////////////////////////////////////////////////////////////////
    
    
    
    
////////////////////////////////////////////////////atop orientation computation/////////////////////////////////////////////////////////////////////////////     
     
     //I use the names of variables from H bonds but has nothing to do with water molecules - just a quick hack
   if(atop_orient=="yes" && type=="density-top")
   {
   // cout<<"Hungarian Horntail"<<endl;
  int *h2onum;
  int *h2oneigh;
  float *h2oxyz;
  float *h2oneighdist;
  float *h2odistmat, *dev_h2odistmat;
  
  //int OHsamples=count_metal;
    int OHsamples=20;
    h2onum = (int *) malloc (sizeof(int)*nElt1*(nstruct/skip));
    h2oneigh = (int *) malloc (sizeof(int)*OHsamples*nElt1*(nstruct/skip));
    h2oneighdist = (float *) malloc (sizeof(float)*OHsamples*nElt1*(nstruct/skip));
    h2oxyz = (float *) malloc (sizeof(float)*7*nElt1*(nstruct/skip));  //Oxyz,H1xyz,H2xyz,Hmidpoint,dipole alignment with surface normal
    h2odistmat = (float *) malloc (sizeof(float)*nElt1*count_metal);
    
  for(int i=0;i<nstruct/skip;i++)
  {
    for(int j=0; j<nElt1; j++)
    {
      h2onum[j+i*nElt1] =0;
      for(int k =0 ;k<OHsamples;k++)
      {
      h2oneigh[j*OHsamples+i*nElt1*OHsamples+k]=0;
      h2oneighdist[j*OHsamples+i*nElt1*OHsamples+k]=0.0;
      }
      for(int k =0 ;k<7;k++)
      {
      h2oxyz[j*7+i*nElt1*7+k]=0.0;
      }
    }
  }   
      
      
cout <<"gpu begins OH distance computation"<<endl;
//Allocate memory in GPU device

    cudaMalloc((void **)&dev_A,sizeof(float)*nElt1*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_B,sizeof(float)*count_metal*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_h2odistmat,sizeof(float)*nElt1*count_metal);
    cudaMalloc((void **)&dev_lattice,sizeof(float)*6*(nstruct/skip));
//Copy data from host to device

    cudaMemcpy(dev_A,A,sizeof(float)*nElt1*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_B,METAL,sizeof(float)*count_metal*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_lattice,lattice,sizeof(int)*6*(nstruct/skip),cudaMemcpyHostToDevice);
//dim3 struct to define elements of the execution configuration

    dim3 dimBlock(32,32,1);
    dim3 dimGrid((nElt1+dimBlock.x-1)/dimBlock.x,(count_metal+dimBlock.y-1)/dimBlock.y,1);
    

   
    
for(int i = 0; i<(nstruct/skip); i++) // goint through every frame
{
  
 // neigh << i<<endl<<endl;
  
   for(int i1=0;i1<nElt1;i1++) // initializing the distance matrix
  {
    for(int j1=0; j1<count_metal; j1++)
    {
      h2odistmat[j1+i1*count_metal]=0.0;
    }
  }
  
  cudaMemcpy(dev_h2odistmat,h2odistmat,sizeof(float)*nElt1*count_metal,cudaMemcpyHostToDevice);
  //All O atoms within bondistOM distance is selected and stored in the above matrix
  
    if(cell_type == "orthorhombic")
    {
          //Hbondmat<<<dimGrid,dimBlock>>>(dev_A,dev_B,dev_h2odistmat,bondistOM,nElt1,count_metal,dev_lattice,i);
	  Hbondmat<<<dimGrid,dimBlock>>>(dev_A,dev_B,dev_h2odistmat,bondistOM,nElt1,count_metal,dev_lattice,i,dirn);
    }
    else if(cell_type == "monoclinic")
    {
          Hbondmatmono<<<dimGrid,dimBlock>>>(dev_A,dev_B,dev_h2odistmat,bondistOM,nElt1,count_metal,dev_lattice,i,dirn);  
	  //Hbondmatmono<<<dimGrid,dimBlock>>>(dev_A,dev_B,dev_h2odistmat,bondistOM,nElt1,count_metal,dev_lattice,i); 
    }

cudaMemcpy(h2odistmat,dev_h2odistmat,sizeof(float)*nElt1*count_metal,cudaMemcpyDeviceToHost);
   
   //cout<<count_metal<<endl;
   
  for(int i1=0;i1<nElt1;i1++) //for each O atoms
  {
    int k=0;
    for(int j1=0; j1<count_metal; j1++) // for each Cu atoms
    {
      if(h2odistmat[j1+i1*count_metal] !=0.0)
      {
	h2onum[i1+i*nElt1]++;
	h2oneigh[i1*OHsamples+i*nElt1*OHsamples+k] = j1;
	h2oneighdist[i1*OHsamples+i*nElt1*OHsamples+k] = h2odistmat[j1+i1*count_metal];
	k++;
      }
    }
  }

}   //frames done

// neigh.close();
    cudaFree(dev_A);
    cudaFree(dev_B);
    cudaFree(dev_h2odistmat);
    cudaFree(dev_lattice);
//sorting needs to be done!!!!

cout<<"Swapping of bonds based on distance "<<endl;

float swapa,swapb;
    for(int i=0;i<(nstruct/skip);i++)
     {
      for(int j=0; j<nElt1; j++)
       {
	 
	for(int swf=0; swf<h2onum[j+i*nElt1]-1;swf++)
	{
	for(int sws=swf+1; sws<h2onum[j+i*nElt1];sws++)
	{
	if(h2oneighdist[j*OHsamples+i*nElt1*OHsamples+swf] > h2oneighdist[j*OHsamples+i*nElt1*OHsamples+sws] && sws !=swf && h2oneighdist[j*OHsamples+i*nElt1*OHsamples+sws] !=0.0 && h2oneighdist[j*OHsamples+i*nElt1*OHsamples+swf] !=0.0) 
	{
	   swapa = h2oneigh[j*OHsamples+i*nElt1*OHsamples+swf]; 
	    h2oneigh[j*OHsamples+i*nElt1*OHsamples+swf]=h2oneigh[j*OHsamples+i*nElt1*OHsamples+sws];
	    h2oneigh[j*OHsamples+i*nElt1*OHsamples+sws]=swapa;
	   swapb = h2oneighdist[j*OHsamples+i*nElt1*OHsamples+swf]; 
	    h2oneighdist[j*OHsamples+i*nElt1*OHsamples+swf]=h2oneighdist[j*OHsamples+i*nElt1*OHsamples+sws];
	    h2oneighdist[j*OHsamples+i*nElt1*OHsamples+sws]=swapb;
	}
	}
	}	
       }
     }	
   


   

cout <<"gpu ends"<<endl;
   
cout<<" done !!!"<<endl;    
    

//h2oneigh contain the IDs of all possible Cu neighbour atoms.


  
  
  float *angdist;
  
  angdist = (float *) malloc (sizeof(float)*180);
  
  for(int i=0; i<=180;i++)
  {
    angdist[i]=0.0;
  }
  
  float val;
  ofstream atopjdist;
  atopjdist.open("atopjdist.dat");
for(int i=0; i<(nstruct/skip); i++) // each frame
{
  //cout<<i<<endl;
 for(int j=0; j<nElt1;j++) // for each oxygen
 {
    if(order=="YZX") {
      val=A[j*3+i*nElt1*3+2];    
    }
    if(order=="XYZ") {
      val=A[j*3+i*nElt1*3+1];      
    }
    if(order=="ZXY") {
      val=A[j*3+i*nElt1*3];      
    }
    //cout<<j<<endl;
   if(h2onum[j+i*nElt1]>0 && exch[j+i*nElt1]==1 && val < maxtop && val > mintop)
   {
     for(int k=0; k<3;k++)
     {
     h2oxyz[j*7+i*nElt1*7+k]= A[j*3+i*nElt1*3+k]; // oxygen atom position
     }
     for(int k=0; k<3;k++)
     {
     h2oxyz[j*7+i*nElt1*7+3+k]= METAL[h2oneigh[j*OHsamples+i*nElt1*OHsamples]*3+i*count_metal*3+k];//Cu position
     //watercount << h2oxyz[j*7+i*nElt1*7+3+k] <<" ";
     }
     //check_together(h2oxyz[j*7+i*nElt1*7],h2oxyz[j*7+i*nElt1*7+1],h2oxyz[j*7+i*nElt1*7+2],\
                    h2oxyz[j*7+i*nElt1*7+3],h2oxyz[j*7+i*nElt1*7+4],h2oxyz[j*7+i*nElt1*7+5],\
                    lattice[i*6],lattice[i*6+1],lattice[i*6+2],lattice[i*6+3],lattice[i*9+4],lattice[i*6+5]);
                    
                    
     float chk1,chk2,chk3;
     float a1,b1,c1,d1,e1,f1;
     float latx,laty,latz,latxy,latxz,latyz;
     //check_together(h2oxyz[j*15+i*nElt1*15],h2oxyz[j*15+i*nElt1*15+1],h2oxyz[j*15+i*nElt1*15+2],\
                    h2oxyz[j*15+i*nElt1*15+3],h2oxyz[j*15+i*nElt1*15+4],h2oxyz[j*15+i*nElt1*15+5],\
                    h2oxyz[j*15+i*nElt1*15+6],h2oxyz[j*15+i*nElt1*15+7],h2oxyz[j*15+i*nElt1*15+8],\
                    lattice[i*6],lattice[i*6+1],lattice[i*6+2],lattice[i*6+3],lattice[i*9+4],lattice[i*6+5]);
     
     a1=h2oxyz[j*7+i*nElt1*7];
     b1=h2oxyz[j*7+i*nElt1*7+1];
     c1=h2oxyz[j*7+i*nElt1*7+2];
     d1=h2oxyz[j*7+i*nElt1*7+3];
     e1=h2oxyz[j*7+i*nElt1*7+4];
     f1=h2oxyz[j*7+i*nElt1*7+5];

     latx=lattice[i*6];
     laty=lattice[i*6+1];
     latz=lattice[i*6+2];
     latxy=lattice[i*6+3];
     latxz=lattice[i*6+4];
     latyz=lattice[i*6+5];
     
chk1= (a1-d1);
chk2= (b1-e1);
chk3= (c1-f1);


if(cell_type == "orthorhombic")
{
if(fabs(chk1) > latx/2.0) { if (chk1 >0) d1=d1+latx; else d1=d1-latx;}
if(fabs(chk2) > laty/2.0) { if (chk2 >0) e1=e1+laty; else e1=e1-laty;}
if(fabs(chk3) > latz/2.0) { if (chk3 >0) f1=f1+latz; else f1=f1-latz;}
}

else if(cell_type == "monoclinic")
{
if(fabs(chk2) > laty/2.0) { if (chk2 >0) {e1=e1+laty; d1=d1+latxy;} else {e1=e1-laty; d1=d1-latxy;}}
chk1= (a1-d1);
if(fabs(chk1) > latx/2.0) { if (chk1 >0) d1=d1+latx; else d1=d1-latx;}
if(fabs(chk3) > latz/2.0) { if (chk3 >0) f1=f1+latz; else f1=f1-latz;}
}
     
     h2oxyz[j*7+i*nElt1*7+3]=d1;
     h2oxyz[j*7+i*nElt1*7+4]=e1;
     h2oxyz[j*7+i*nElt1*7+5]=f1;
                    
                    
                    
     //////////////////////////////angle between dipole and z axis (specific - be careful if using other directions//////////////////////////////////////////////////
     //h2oxyz[j*7+i*nElt1*7+6]=angle(h2oxyz[j*7+i*nElt1*7+3],h2oxyz[j*7+i*nElt1*7+4],h2oxyz[j*7+i*nElt1*7+5],\
                                   h2oxyz[j*7+i*nElt1*7+3],h2oxyz[j*7+i*nElt1*7+4],h2oxyz[j*7+i*nElt1*7+5]+5.0,\
                                   h2oxyz[j*7+i*nElt1*7+3],h2oxyz[j*7+i*nElt1*7+4],h2oxyz[j*7+i*nElt1*7+5],\
                                   h2oxyz[j*7+i*nElt1*7],h2oxyz[j*7+i*nElt1*7+1],h2oxyz[j*7+i*nElt1*7+2]);
       //h2oxyz[j*7+i*nElt1*7+6]=angle(0,0,0,0,0,1,\
                                   h2oxyz[j*7+i*nElt1*7+3],h2oxyz[j*7+i*nElt1*7+4],h2oxyz[j*7+i*nElt1*7+5],\
                                   h2oxyz[j*7+i*nElt1*7],h2oxyz[j*7+i*nElt1*7+1],h2oxyz[j*7+i*nElt1*7+2]);
       h2oxyz[j*7+i*nElt1*7+6]=angle(h2oxyz[j*7+i*nElt1*7+3],h2oxyz[j*7+i*nElt1*7+4],h2oxyz[j*7+i*nElt1*7+5],\
                                    h2oxyz[j*7+i*nElt1*7+3],h2oxyz[j*7+i*nElt1*7+4],h2oxyz[j*7+i*nElt1*7+5]+1,\
                                   h2oxyz[j*7+i*nElt1*7+3],h2oxyz[j*7+i*nElt1*7+4],h2oxyz[j*7+i*nElt1*7+5],\
                                   h2oxyz[j*7+i*nElt1*7],h2oxyz[j*7+i*nElt1*7+1],h2oxyz[j*7+i*nElt1*7+2]);
     //watercount << h2oxyz[j*7+i*nElt1*7+3] <<" "<< h2oxyz[j*7+i*nElt1*7+3+1] <<" "<< h2oxyz[j*7+i*nElt1*7+3+2] <<" ";
     //watercount << h2oxyz[j*7+i*nElt1*7+6]*(180.0/3.141);
     angdist[int(floor(h2oxyz[j*7+i*nElt1*7+6]*(180.0/3.14159)))]++;
    // watercount<<endl;
     atopjdist << j<<" "<<h2oxyz[j*7+i*nElt1*7+6]*(180.0/3.14159)<<endl;
   }
   
  }
}
  //watercount.close();
  atopjdist.close();
float totalatopdist=0.0;
  for(int i=0; i<=180;i++)
  {
    totalatopdist +=angdist[i];
  }
  ofstream angldist;
  angldist.open("atop-orient.dat");
  angldist <<"T num"<<endl;
   for(int i=0; i<=180;i++)
  {
    //angldist << i<<" "<<angdist[i]/(totalatopdist)<<endl;
    angldist << i<<" "<<angdist[i]/float(nstruct/skip)<<endl;
  }
  angldist.close();
  free(angdist);
  free(h2onum);
  free(h2oneigh);
  free(h2oneighdist);
  free(h2oxyz);
  free(h2odistmat);     
     
     
   }
    
////////////////////////////////////////////////////end of atop orientation calculation /////////////////////////////////////////////////////////////////////////////////
    
    
    
    
    
    
       
    
    
    
    
    
    
    
    if(filter_density == "yes")
    {
    cudaMalloc((void **)&dev_densityf,sizeof(int)*xsplit*ysplit*zsplit);
    cudaMemcpy(dev_densityf,densityf,sizeof(int)*xsplit*ysplit*zsplit,cudaMemcpyHostToDevice);
    dim3 dimBlockee(10,10,10);
    dim3 dimGridee((xsplit+dimBlockee.x-1)/dimBlockee.x,(ysplit+dimBlockee.y-1)/dimBlockee.y,(zsplit+dimBlockee.z-1)/dimBlockee.y);
    if(order=="XYZ") {
    density_filter<<<dimGridee,dimBlockee>>>(dev_density, dev_densityf,ysplit,zsplit,xsplit);
    }
    else if(order=="ZXY") {
    density_filter<<<dimGridee,dimBlockee>>>(dev_density, dev_densityf,xsplit,ysplit,zsplit);
    }
    else if(order=="YZX") {
    density_filter<<<dimGridee,dimBlockee>>>(dev_density, dev_densityf,zsplit,xsplit,ysplit);
    }
    cudaMemcpy(density,dev_densityf,sizeof(int)*xsplit*ysplit*zsplit,cudaMemcpyDeviceToHost);
    }
    else
    {
      cudaMemcpy(density,dev_density,sizeof(int)*xsplit*ysplit*zsplit,cudaMemcpyDeviceToHost);
    }
FILE *densityfile=fopen("out","wt");



 if(plot =="int" || plot == "bulk")
  {
    int diffdirn;
if(vacuum.at(0)=='z') diffdirn=2;
if(vacuum.at(0)=='y') diffdirn=1;
if(vacuum.at(0)=='x') diffdirn=0;

  float minAz,maxAz,midAz;
  minAz=A_int[diffdirn];
  maxAz=A_int[diffdirn];
  for(int j=0; j<nA_int;j++)
    {
      if(A_int[j*3+diffdirn] > maxAz) maxAz= A_int[j*3+diffdirn];
      if(A_int[j*3+diffdirn] < minAz) minAz= A_int[j*3+diffdirn];
    }
    midAz=(maxAz+minAz)/2.0;
  ofstream diffplot;
  ofstream diffplot1;
  diffplot.open("diffplot");
  diffplot1.open("diffplot1");
  if (type=="density" && scope=="bulk-interface" && plot=="int")
 {
   for(int i=0; i<nstruct/skip;i++)
   {
    for(int j=0; j<nA_int;j++)
    {
      if(A_int[j*3+diffdirn] < midAz)
      {
    diffplot <<A_int[j*3+i*nA_int*3]<<" "<<A_int[j*3+i*nA_int*3+1]<<" "<<A_int[j*3+i*nA_int*3+2]<<endl;
      }
      else 
      {
    diffplot1 <<A_int[j*3+i*nA_int*3]<<" "<<A_int[j*3+i*nA_int*3+1]<<" "<<A_int[j*3+i*nA_int*3+2]<<endl;
      }
    }
   }
 }
 if (type=="density" && scope=="bulk-interface" && plot=="bulk")
 {
   for(int i=0; i<nstruct/skip;i++)
   {
    for(int j=0; j<nA_bulk;j++)
    {
      if(A_bulk[j*3+diffdirn] < midAz)
      {
    diffplot <<A_bulk[j*3+i*nA_bulk*3]<<" "<<A_bulk[j*3+i*nA_bulk*3+1]<<" "<<A_bulk[j*3+i*nA_bulk*3+2]<<endl;
      }
      else 
      {
    diffplot1 <<A_bulk[j*3+i*nA_bulk*3]<<" "<<A_bulk[j*3+i*nA_bulk*3+1]<<" "<<A_bulk[j*3+i*nA_bulk*3+2]<<endl;
      }
    }
   }
 }
 diffplot.close();
 diffplot1.close();  
  }




if(order=="XYZ") {

for(int l=0;l<ysplit;l++)
{
  //fprintf(densityfile,"%d %d\n",zsplit,xsplit);
  fprintf(densityfile,"#%d \n",l);
for(int m=0;m<zsplit;m++)
{
for(int k=0;k<xsplit;k++)
{
    if(type=="density-top") 
  {
fprintf(densityfile,"%f ",float(density[m+l*zsplit+k*ysplit*zsplit])/float(maxbondist-minbondist));
  }
  else
  {
    fprintf(densityfile,"%f ",float(density[m+l*zsplit+k*ysplit*zsplit])/float(nElt1));
  }
}
fprintf(densityfile,"\n");
}
//fprintf(densityfile,"5\n");
//fprintf(densityfile,"1 10.\n");
//fprintf(densityfile,"2 15.\n");
//fprintf(densityfile,"3 20.\n");
//fprintf(densityfile,"4 25.\n");
//fprintf(densityfile,"5 30.\n");
//fprintf(densityfile,"6 \n"); 
}
}
if(order=="ZXY") {

for(int k=0;k<xsplit;k++)
{
  fprintf(densityfile,"#%d \n",k);
for(int l=0;l<ysplit;l++)
{
for(int m=0;m<zsplit;m++)
{
      if(type=="density-top") 
  {
fprintf(densityfile,"%f ",float(density[m+l*zsplit+k*ysplit*zsplit])/float(maxbondist-minbondist));
  }
    else
  {
    fprintf(densityfile,"%f ",float(density[m+l*zsplit+k*ysplit*zsplit])/float(nElt1));
  }
}
fprintf(densityfile,"\n");
}
//fprintf(densityfile,"5\n");
//fprintf(densityfile,"1 10.\n");
//fprintf(densityfile,"2 15.\n");
//fprintf(densityfile,"3 20.\n");
//fprintf(densityfile,"4 25.\n");
//fprintf(densityfile,"5 30.\n");
//fprintf(densityfile,"6 \n");
}
}
if(order=="YZX") {

for(int m=0;m<zsplit;m++)
{
  fprintf(densityfile,"#%d \n",m);
for(int k=0;k<xsplit;k++)
{
for(int l=0;l<ysplit;l++)
{
if(type=="density-top") 
  {
fprintf(densityfile,"%f ",float(density[m+l*zsplit+k*ysplit*zsplit])/float(maxbondist-minbondist));
  }
else
  {
fprintf(densityfile,"%f ",float(density[m+l*zsplit+k*ysplit*zsplit])/float(nElt1));
  }
}
fprintf(densityfile,"\n");
}
//fprintf(densityfile,"5\n");
//fprintf(densityfile,"1 10.\n");
//fprintf(densityfile,"2 15.\n");
//fprintf(densityfile,"3 20.\n");
//fprintf(densityfile,"4 25.\n");
//fprintf(densityfile,"5 30.\n");
//fprintf(densityfile,"6 \n");
}  
}
fclose(densityfile);

cudaFree(dev_density);
cudaFree(dev_densityf);
if (type=="density" || type=="densityvelocity");
{
cudaFree(dev_A);
cudaFree(dev_xtick);
cudaFree(dev_ytick);
cudaFree(dev_ztick);
}

if(write_cube == "yes")
{
  FILE *densityfile=fopen("cube","wt");
  if(order=="XYZ") {

for(int l=0;l<ysplit;l++)
{
for(int m=0;m<zsplit;m++)
{
for(int k=0;k<xsplit;k++)
{
fprintf(densityfile,"%d ",density[m+l*zsplit+k*ysplit*zsplit]);
}
}
}
  }
else if(order=="ZXY") {

for(int k=0;k<xsplit;k++)
{
for(int l=0;l<ysplit;l++)
{
for(int m=0;m<zsplit;m++)
{
fprintf(densityfile,"%f ",float(density[m+l*zsplit+k*ysplit*zsplit]));
if (m % 6 == 5) fprintf(densityfile,"\n");
}
fprintf(densityfile,"\n");
}
}
}
else if(order=="YZX") {

for(int k=0;k<xsplit;k++)
{
for(int l=0;l<ysplit;l++)
{
for(int m=0;m<zsplit;m++)
{
fprintf(densityfile,"%f ",float(density[m+l*zsplit+k*ysplit*zsplit]));
if (m % 6 == 5) fprintf(densityfile,"\n");
}
fprintf(densityfile,"\n");
}
}
}

fclose(densityfile);
}

if(type=="density-top" || (type=="density" && scope=="bulk-interface"))
{
  system("paste out out out | tail -n +2 > ext-out ; cat ext-out ext-out ext-out > extended.out ");
}

}
