//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"

void density_side_avg_pos()
{


     for(int i=0; i< nstruct/skip;i=i++)
     {
      for(int j=0; j< nElt1; j++)
      {
       if(i>0)
       {
        if(fabs(A[j*3+i*nElt1*3]-A[j*3]) > xvec/2)
        {
         A[j*3+i*nElt1*3]=-A[j*3+i*nElt1*3]+xvec;
        }
        if(fabs(A[j*3+i*nElt1*3+1]-A[j*3+1]) > yvec/2)
        {
         A[j*3+i*nElt1*3+1]=-A[j*3+i*nElt1*3+1]+yvec;
        }
        if(fabs(A[j*3+i*nElt1*3+2]-A[j*3+2]) > zvec/2)
        {
         A[j*3+i*nElt1*3+2]=-A[j*3+i*nElt1*3+2]+zvec;
        }
       }
       Aavg[j*3]+=A[j*3+i*nElt1*3];
       Aavg[j*3+1]+=A[j*3+i*nElt1*3+1];
       Aavg[j*3+2]+=A[j*3+i*nElt1*3+2];
      }
     }
      for(int j=0; j< nElt1; j++)
      {     
       Aavg[j*3]=Aavg[j*3]/(nstruct/skip);
       Aavg[j*3+1]=Aavg[j*3+1]/(nstruct/skip);
       Aavg[j*3+2]=Aavg[j*3+2]/(nstruct/skip); 
     //  printf("%f %f %f \n",Aavg[j*3],Aavg[j*3+1],Aavg[j*3+2]);
      }
     nstruct=1;
    
    FILE *testout=fopen("teter.dat","wt");
      for(int j=0; j< nElt1; j++)
      {
	if(order=="YZX") {
	fprintf(testout,"%s \n",nElt1);
	fprintf(testout,"%d \n", zsplit);
       fprintf(testout,"%d %f %d %f %f \n",symbol,size,color,Aavg[j*3+1]/yvec,Aavg[j*3]/xvec);       
	}
	if(order=="XYZ") {
	fprintf(testout,"%s \n",nElt1);
	fprintf(testout,"%d \n", ysplit);  
       fprintf(testout,"%d %f %d %f %f \n",symbol,size,color,Aavg[j*3]/xvec,Aavg[j*3+2]/zvec);
 //printf("%d %f %d %f %f %f\n",symbol,size,color,Aavg[j*3],Aavg[j*3+1],Aavg[j*3+2]);      
	}
	if(order=="ZXY") {
  	fprintf(testout,"%s \n",nElt1);
	fprintf(testout,"%d \n", xsplit);
       fprintf(testout,"%d %f %d %f %f \n",symbol,size,color,Aavg[j*3+2]/zvec,Aavg[j*3+1]/yvec);
      // printf("%d %f %d %f %f %f\n",symbol,size,color,Aavg[j*3],Aavg[j*3+1],Aavg[j*3+2]); 
	}
      }    
    fclose(testout);

}