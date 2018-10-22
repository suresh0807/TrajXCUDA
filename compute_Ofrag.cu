//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################


//both Elt1 and Elt2 are O atoms. The problem here is that O-O distances are only used to find H bonds,
//which can be very wrong. Check compute_hbonds instead.

#include "cudatools.cuh"
int *Ofrag;
int *Ofragsum;
int *dev_Ofrag;



void compute_Ofrag(void)
{
 Ofrag= (int *) malloc (sizeof(int)*nElt1*(nstruct/skip));
// Initialize distance matrix and histogram matrix
 for(int l=0; l<nstruct/skip; l++)
    {
      for(int i=0; i< nElt1; i++)
        {
          Ofrag[i+l*nElt1]=0;
        }
      }
 Ofragsum= (int *) malloc (sizeof(int)*4*(nstruct/skip));   
 for(int l=0; l<nstruct/skip; l++)
    {
      for(int i=0; i< 4; i++)
        {
          Ofragsum[i+l*4]=0;
        }
      }
      

//Allocate memory in GPU device

    cudaMalloc((void **)&dev_A,sizeof(float)*nElt1*nstruct/skip*3);
    cudaMalloc((void **)&dev_B,sizeof(float)*nElt2*nstruct/skip*3);
    cudaMalloc((void **)&dev_Ofrag,sizeof(int)*nElt1*nstruct);
//Copy data from host to device

    cudaMemcpy(dev_A,A,sizeof(float)*nElt1*nstruct/skip*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_B,B,sizeof(float)*nElt2*nstruct/skip*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_Ofrag,Ofrag,sizeof(int)*nElt1*nstruct/skip,cudaMemcpyHostToDevice);
//dim3 struct to define elements of the execution configuration


    dim3 dimBlock(10,10,10);
    dim3 dimGrid((nElt1+dimBlock.x-1)/dimBlock.x,(nElt2+dimBlock.y-1)/dimBlock.y,((nstruct/skip)+dimBlock.z-1)/dimBlock.z);

//Cuda kernal execution for distance matrix with CUDA timing API commands

    Ofragmat<<<dimGrid,dimBlock>>>(dev_A,dev_B,dev_Ofrag,bondist,nElt1,nElt2,nstruct/skip,xvec,yvec,zvec);

    cudaMemcpy(Ofrag,dev_Ofrag,sizeof(int)*nElt1*nstruct/skip,cudaMemcpyDeviceToHost);
    
    cudaFree(dev_A);
    cudaFree(dev_B);
    
 for(int l=0; l<nstruct/skip; l++)
    {
       for(int i=0; i< nElt1; i++)
         {
          if(Ofrag[i+l*nElt1] ==2)
	  {
	    Ofragsum[l*4]+=1;
	  }
	  else if(Ofrag[i+l*nElt1] ==1)
	  {
	    Ofragsum[(l*4)+1]+=1;
	  }
	  else if(Ofrag[i+l*nElt1] ==0)
	  {
	    Ofragsum[(l*4)+2]+=1;
	  }
	  else 
	  {
	    Ofragsum[(l*4)+3]+=1;
	  }	   
	}
     }
    
    
    

//    FILE *Ofragplotchk=fopen("O-fragments-chk.data","wt");
// for(int l=0; l<nstruct; l++)
//    {
//      fprintf(Ofragplotchk,"%d\n",l);
//      for(int i=0; i< nElt1; i++)
//        {
//          fprintf(Ofragplotchk,"%d %d\n",i,Ofrag[i+l*nElt1]);
//        }
//      }
        
//    fclose(Ofragplotchk);
    free(Ofrag);
    
    FILE *Ofragplot=fopen("O-fragments.data","wt");
    fprintf(Ofragplot,"#frame h2o oh o h3o\n");
 for(int l=0; l<nstruct/skip; l++)
    {
          fprintf(Ofragplot,"%d %d %d %d %d %d %d\n",l,Ofragsum[l*4],Ofragsum[1+l*4],Ofragsum[2+l*4],Ofragsum[3+l*4],Ofragsum[l*4]+Ofragsum[1+l*4]+Ofragsum[2+l*4]+Ofragsum[3+l*4],(Ofragsum[l*4]+Ofragsum[1+l*4]+Ofragsum[2+l*4]+Ofragsum[3+l*4])*2);
       }
        //frame h2o oh o h3o
    fclose(Ofragplot);

  free(Ofragsum);
}
