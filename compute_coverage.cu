//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################




#include "cudatools.cuh"
int *Mcover, *Scover;
int *Mcoversum, *Scoversum;
int *dev_cover;


void compute_coverage(void)
{
 Mcover= (int *) malloc (sizeof(int)*nElt1*(nstruct/skip));
 Scover= (int *) malloc (sizeof(int)*nElt2*(nstruct/skip));
// Initialize distance matrix and histogram matrix
 for(int l=0; l<nstruct/skip; l++)
    {
      for(int i=0; i< nElt1; i++)
        {
          Mcover[i+l*nElt1]=0;
        }
      for(int i=0; i< nElt2; i++)
        {
          Scover[i+l*nElt2]=0;
        }        
      }

 Mcoversum= (int *) malloc (sizeof(int)*(nstruct/skip));   
 Scoversum= (int *) malloc (sizeof(int)*(nstruct/skip));
 for(int l=0; l<nstruct/skip; l++)
    {
          Mcoversum[l]=0;
	  Scoversum[l]=0;
      }
//Allocate memory in GPU device

    cudaMalloc((void **)&dev_A,sizeof(float)*nElt1*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_B,sizeof(float)*nElt1*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_cover,sizeof(int)*nElt1*(nstruct/skip));


//Copy data from host to device

    cudaMemcpy(dev_A,A,sizeof(float)*nElt1*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dev_B,A,sizeof(float)*nElt1*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    

//dim3 struct to define elements of the execution configuration


    dim3 dimBlock(10,10,10);
    dim3 dimGrid((nElt1+dimBlock.x-1)/dimBlock.x,(nElt1+dimBlock.y-1)/dimBlock.y,((nstruct/skip)+dimBlock.z-1)/dimBlock.z);

//Cuda kernal execution for distance matrix with CUDA timing API commands
 
    //For surface metal atoms
    cudaMemcpy(dev_cover,Mcover,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyHostToDevice);
    if (cell_type == "orthorhombic") 
    {
    covermat<<<dimGrid,dimBlock>>>(dev_A,dev_B,dev_cover,Mbondist,nElt1,nElt1,(nstruct/skip),xvec,yvec,zvec);
    }
    else if (cell_type == "monoclinic")
    {
     covermatmono<<<dimGrid,dimBlock>>>(dev_A,dev_B,dev_cover,Mbondist,nElt1,nElt1,(nstruct/skip),xvec,yvec,zvec,yxvec);
    }
    cudaMemcpy(Mcover,dev_cover,sizeof(int)*nElt1*(nstruct/skip),cudaMemcpyDeviceToHost);
    //For substrate atoms
    cudaFree(dev_B);
    cudaFree(dev_cover);
    
    
    
    cudaMalloc((void **)&dev_B,sizeof(float)*nElt2*(nstruct/skip)*3);
    cudaMalloc((void **)&dev_cover,sizeof(int)*nElt2*(nstruct/skip));
    cudaMemcpy(dev_B,B,sizeof(float)*nElt2*(nstruct/skip)*3,cudaMemcpyHostToDevice);
    
    dim3 dimBlocks(10,10,10);
    dim3 dimGrids((nElt1+dimBlock.x-1)/dimBlock.x,(nElt2+dimBlock.y-1)/dimBlock.y,((nstruct/skip)+dimBlock.z-1)/dimBlock.z);
    
    cudaMemcpy(dev_cover,Scover,sizeof(int)*nElt2*(nstruct/skip),cudaMemcpyHostToDevice);
    if (cell_type == "orthorhombic") 
    {
    covermat<<<dimGrids,dimBlocks>>>(dev_B,dev_A,dev_cover,bondist,nElt2,nElt1,(nstruct/skip),xvec,yvec,zvec);
    }
    else if (cell_type == "monoclinic")
    {
     covermatmono<<<dimGrids,dimBlocks>>>(dev_B,dev_A,dev_cover,bondist,nElt2,nElt1,(nstruct/skip),xvec,yvec,zvec,yxvec);
    }
    cudaMemcpy(Scover,dev_cover,sizeof(int)*nElt2*(nstruct/skip),cudaMemcpyDeviceToHost);

    cudaFree(dev_A);
    cudaFree(dev_B);
    cudaFree(dev_cover);
    
    
     for(int l=0; l<(nstruct/skip); l++)
    {
       for(int i=0; i< nElt1; i++)
         {
	   //cout <<l<<" "<<i<<" "<<Mcover[i+l*nElt1]<<endl;
          if(Mcover[i+l*nElt1] < maxneigh && Mcover[i+l*nElt1] > minneigh)
	  {
	    Mcoversum[l]+=1;
	  }
	}
    }
    
    for(int l=0; l<(nstruct/skip); l++)
    {
       for(int i=0; i< nElt2; i++)
         {
	  // cout <<l<<" "<<i<<" "<<Scover[i+l*nElt2]<<endl;
          if( Scover[i+l*nElt2] > 0)
	  {
	    Scoversum[l]+=1;
	  }
	}
    }
    
    int *sorter = (int *) malloc (sizeof(int)*nElt1*(nstruct/skip));

    int surf_atom_counter=0;
    float int_max, int_min;
    
     for(int l=0; l<(nstruct/skip); l++)
    {
       for(int i=0; i< nElt1; i++)
         {
          if(Mcover[i+l*nElt1] < maxneigh && Mcover[i+l*nElt1] > minneigh){sorter[i+l*nElt1] = 0;if(l==0){surf_atom_counter++;}}
          else {sorter[i+l*nElt1]=1;}
          //cout << l<<" "<<i<<" "<< sorter[i+l*nElt1] <<" "<< A[i*3+l*nElt1*3+2]<<endl;        
	}
    }    
    cout <<surf_atom_counter<<endl;
    float *int_maxmin = (float *) malloc (sizeof(float)*surf_atom_counter*3);
    
       int counta=0;
       for(int i=0; i< nElt1; i++)
       {
          if(sorter[i] == 0){int_maxmin[counta] = A[i*3];
	                     int_maxmin[counta+1] = A[i*3+1];
			     int_maxmin[counta+2] = A[i*3+2];
			     counta+=3;
	    // cout <<A[i*3+2]<<endl;
	  }       
       }    
    
    int_max = maximum(int_maxmin,surf_atom_counter,2);
    int_min = minimum(int_maxmin,surf_atom_counter,2);
    
    cout <<"Minimum Z position "<<int_min<<endl;
    cout <<"Maximum Z position "<<int_max<<endl;
    
    
    
// printing the rdf data to be visualized
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
    free(Mcover);
    free(Scover);
    ofstream coverplot;
    coverplot.open("percentage-coverage.data");
    //FILE *coverplot=fopen("percentage-coverage.data","wt");
 for(int l=0; l<(nstruct/skip); l++)
    {
          coverplot << l <<" "<< Mcoversum[l] <<" "<< Scoversum[l] <<" "<< (float(Scoversum[l])/float(Mcoversum[l]))*100.0 <<endl;
    }
    coverplot.close();    
    //fclose(coverplot);
  free(Mcoversum);
  free(Scoversum);
}
