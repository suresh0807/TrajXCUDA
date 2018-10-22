//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################


#include "cudatools.cuh"

void density_top_avg_pos()
{

cout <<" Copper average positions - computing"<<endl;

  if(cell_type=="orthorhombic")
  {
     for(int i=0; i< nstruct/skip;i++)
     {
      for(int j=0; j< nElt1; j++)
      {
       if(i>0)
       {
        if(fabs(A[j*3+i*nElt1*3]-A[j*3]) > xvec/2)
        {
	 if(A[j*3+i*nElt1*3]-A[j*3] > 0){ A[j*3+i*nElt1*3]=A[j*3+i*nElt1*3]-xvec;}
	 else if (A[j*3+i*nElt1*3]-A[j*3] < 0){ A[j*3+i*nElt1*3]=A[j*3+i*nElt1*3]+xvec;}
        }
        if(fabs(A[j*3+i*nElt1*3+1]-A[j*3+1]) > yvec/2)
        {
         if(A[j*3+i*nElt1*3+1]-A[j*3+1] > 0){A[j*3+i*nElt1*3+1]=A[j*3+i*nElt1*3+1]-yvec;}
         else if(A[j*3+i*nElt1*3+1]-A[j*3+1] < 0){A[j*3+i*nElt1*3+1]=A[j*3+i*nElt1*3+1]+yvec;}
        }
        if(fabs(A[j*3+i*nElt1*3+2]-A[j*3+2]) > zvec/2)
        {
         if(A[j*3+i*nElt1*3+2]-A[j*3+2] > 0){ A[j*3+i*nElt1*3+2]=A[j*3+i*nElt1*3+2]-zvec;}
         else if(A[j*3+i*nElt1*3+2]-A[j*3+2] < 0){ A[j*3+i*nElt1*3+2]=A[j*3+i*nElt1*3+2]+zvec;}
        }
       }
       Aavg[j*3]+=A[j*3+i*nElt1*3];
       Aavg[j*3+1]+=A[j*3+i*nElt1*3+1];
       Aavg[j*3+2]+=A[j*3+i*nElt1*3+2];
      }
     }
  }
  else if(cell_type=="monoclinic")
  {
    cout<<xvec<<" "<<yvec<<" "<<zvec<<" "<<yxvec<<endl;
     
    for(int i=0; i< nstruct/skip;i++)
     {
      for(int j=0; j< nElt1; j++)
      {
       if(i>0)
       {
	 if(fabs(A[j*3+i*nElt1*3+1]-A[j*3+1]) > yvec/2)
        {
         if(A[j*3+i*nElt1*3+1]-A[j*3+1] > 0){A[j*3+i*nElt1*3+1]=A[j*3+i*nElt1*3+1]-yvec;A[j*3+i*nElt1*3]=A[j*3+i*nElt1*3]-yxvec;}
         else if(A[j*3+i*nElt1*3+1]-A[j*3+1] < 0){A[j*3+i*nElt1*3+1]=A[j*3+i*nElt1*3+1]+yvec;A[j*3+i*nElt1*3]=A[j*3+i*nElt1*3]+yxvec;}
        }
        if(fabs(A[j*3+i*nElt1*3]-A[j*3]) > xvec/2)
        {
	 if(A[j*3+i*nElt1*3]-A[j*3] > 0){ A[j*3+i*nElt1*3]=A[j*3+i*nElt1*3]-xvec;}
	 else if (A[j*3+i*nElt1*3]-A[j*3] < 0){ A[j*3+i*nElt1*3]=A[j*3+i*nElt1*3]+xvec;}
        }
        if(fabs(A[j*3+i*nElt1*3+2]-A[j*3+2]) > zvec/2)
        {
         if(A[j*3+i*nElt1*3+2]-A[j*3+2] > 0){ A[j*3+i*nElt1*3+2]=A[j*3+i*nElt1*3+2]-zvec;}
         else if(A[j*3+i*nElt1*3+2]-A[j*3+2] < 0){ A[j*3+i*nElt1*3+2]=A[j*3+i*nElt1*3+2]+zvec;}
        }
        
        
        if(fabs(A[j*3+i*nElt1*3+1]-A[j*3+1]) > yvec/6)
        {
         if(A[j*3+i*nElt1*3+1]-A[j*3+1] > 0){A[j*3+i*nElt1*3+1]=A[j*3+i*nElt1*3+1]-(yvec/3);A[j*3+i*nElt1*3]=A[j*3+i*nElt1*3]-(yxvec/3);}
         else if(A[j*3+i*nElt1*3+1]-A[j*3+1] < 0){A[j*3+i*nElt1*3+1]=A[j*3+i*nElt1*3+1]+(yvec/3);A[j*3+i*nElt1*3]=A[j*3+i*nElt1*3]+(yxvec/3);}
        }
        if(fabs(A[j*3+i*nElt1*3]-A[j*3]) > xvec/6)
        {
	 if(A[j*3+i*nElt1*3]-A[j*3] > 0){ A[j*3+i*nElt1*3]=A[j*3+i*nElt1*3]-(xvec/3);}
	 else if (A[j*3+i*nElt1*3]-A[j*3] < 0){ A[j*3+i*nElt1*3]=A[j*3+i*nElt1*3]+(xvec/3);}
        }
        if(fabs(A[j*3+i*nElt1*3+2]-A[j*3+2]) > zvec/6)
        {
         if(A[j*3+i*nElt1*3+2]-A[j*3+2] > 0){ A[j*3+i*nElt1*3+2]=A[j*3+i*nElt1*3+2]-(zvec/3);}
         else if(A[j*3+i*nElt1*3+2]-A[j*3+2] < 0){ A[j*3+i*nElt1*3+2]=A[j*3+i*nElt1*3+2]+(zvec/3);}
        }
       }
       Aavg[j*3]+=A[j*3+i*nElt1*3];
       Aavg[j*3+1]+=A[j*3+i*nElt1*3+1];
       Aavg[j*3+2]+=A[j*3+i*nElt1*3+2];
      }
     }
     
    
  }  
      for(int j=0; j< nElt1; j++)
      {     
       Aavg[j*3]=Aavg[j*3]/(nstruct/skip);
       Aavg[j*3+1]=Aavg[j*3+1]/(nstruct/skip);
       Aavg[j*3+2]=Aavg[j*3+2]/(nstruct/skip); 
     //  printf("%f %f %f \n",Aavg[j*3],Aavg[j*3+1],Aavg[j*3+2]);
      }
/*
      for(int j=0; j< nElt1; j++)
      {     
       Aavg[j*3]=A[j*3];
       Aavg[j*3+1]=A[j*3+1];
       Aavg[j*3+2]=A[j*3+2]; 
     //  printf("%f %f %f \n",Aavg[j*3],Aavg[j*3+1],Aavg[j*3+2]);
      }
      */
   /*    for(int j=0; j< nElt1; j++)
      {     
	for(int k=0; k<3;k++)
	{
          if(Aavg[j*3+k] <0) {Aavg[j*3+k]=1-fabs(Aavg[j*3+k]);}
	}
      }*/
cout <<"Writing the average positions"<<endl;  
FILE *test=fopen("avg-pos.dat","wt");
    FILE *testout=fopen("teter.dat","wt");
      /*for(int j=0; j< nElt1; j++)
      {
	if(order=="YZX") {
	if(Aavg[j*3+2] < maxtop && Aavg[j*3+2] > mintop)
	{cout <<"Alles gut "<< order<<endl;
	cout<<nElt1<<endl;
	cout<<zsplit<<endl;
       printf("%d %f %d %f %f \n",symbol,size,color,Aavg[j*3+1]/yvec,Aavg[j*3]/xvec);
	fprintf(testout,"%s \n",nElt1);
	fprintf(testout,"%d \n", zsplit);
       fprintf(testout,"%d %f %d %f %f \n",symbol,size,color,Aavg[j*3+1]/yvec,Aavg[j*3]/xvec);  
	}cout <<"Alles gut "<< order<<endl;
	}
	else if(order=="XYZ") {
	if(Aavg[j*3+1] < maxtop && Aavg[j*3+1] > mintop)
	{
	fprintf(testout,"%s \n",nElt1);
	fprintf(testout,"%d \n", ysplit);  
       fprintf(testout,"%d %f %d %f %f \n",symbol,size,color,Aavg[j*3]/xvec,Aavg[j*3+2]/zvec);
 //printf("%d %f %d %f %f %f\n",symbol,size,color,Aavg[j*3],Aavg[j*3+1],Aavg[j*3+2]);      
	}
	}
	else if(order=="ZXY") {
	if(Aavg[j*3] < maxtop && Aavg[j*3] > mintop)
	{
  	fprintf(testout,"%s \n",nElt1);
	fprintf(testout,"%d \n", xsplit);
       fprintf(testout,"%d %f %d %f %f \n",symbol,size,color,Aavg[j*3+2]/zvec,Aavg[j*3+1]/yvec);
      // printf("%d %f %d %f %f %f\n",symbol,size,color,Aavg[j*3],Aavg[j*3+1],Aavg[j*3+2]); 
	}
	}
      }*/
      	//fprintf(testout,"%d \n",nElt1);
        //fprintf(testout,"%d \n", ysplit);
        //for(int j=0; j< nElt1; j++)
     // {
	//fprintf(test,"%f %f %f \n",Aavg[j*3],Aavg[j*3+1],Aavg[j*3+2]);
      //}
       // fclose(test);
        
      if(order=="YZX") { if(cell_type=="monoclinic") xvec=xvec+yxvec;
      for(int j=0; j< nElt1; j++)
      {
	if(Aavg[j*3+2] <= maxtop && Aavg[j*3+2] >= mintop)
	{  float app = Aavg[j*3+1]*(float(ysplit)/yvec);
	   float bpp = Aavg[j*3]*(float(xsplit)/xvec);
	   //if (app <0) {app = app+1;}
	   //if (bpp <0) {bpp = bpp+1;}
	   //if (app >1) {app = app-1;}
	   //if (bpp >1) {bpp = bpp-1;}
	 // fprintf(testout,"%d %f %d %f %f \n",symbol,size,color,app,bpp);
	 fprintf(testout,"%f %f \n",(app),(bpp));
	  fprintf(test,"%f %f %f \n",Aavg[j*3],Aavg[j*3+1],Aavg[j*3+2]);
	}
      }
      }
      else if(order=="XYZ") {
	for(int j=0; j< nElt1; j++)
      {
	if(Aavg[j*3+1] < maxtop && Aavg[j*3+1] > mintop)
	{  float app = Aavg[j*3]/xvec;
	   float bpp = Aavg[j*3+2]/zvec;
	   if (app <0) {app = app+1;}
	   if (bpp <0) {bpp = bpp+1;}
	   if (app >1) {app = app-1;}
	   if (bpp >1) {bpp = bpp-1;}
	  //fprintf(testout,"%d %f %d %f %f \n",symbol,size,color,app,bpp);
	  fprintf(testout,"%f %f \n",(app*xsplit)+xsplit,(bpp*zsplit)+zsplit);
	  fprintf(test,"%f %f %f \n",Aavg[j*3],Aavg[j*3+1],Aavg[j*3+2]);
	}
      }
      }
           else if(order=="ZXY") {
	for(int j=0; j< nElt1; j++)
      {
	if(Aavg[j*3] < maxtop && Aavg[j*3] > mintop)
	{  float app = Aavg[j*3+2]/zvec;
	   float bpp = Aavg[j*3+1]/yvec;
	   if (app <0) {app = app+1;}
	   if (bpp <0) {bpp = bpp+1;}
	   if (app >1) {app = app-1;}
	   if (bpp >1) {bpp = bpp-1;}
	  //fprintf(testout,"%d %f %d %f %f \n",symbol,size,color,app,bpp);
	  fprintf(testout,"%f %f \n",(app*zsplit)+zsplit,(bpp*ysplit)+ysplit);
	  fprintf(test,"%f %f %f \n",Aavg[j*3],Aavg[j*3+1],Aavg[j*3+2]);
	}
      }
      }
    fclose(testout);
    fclose(test);

    


    
    
        
cout <<"Alles gut"<<endl;
}