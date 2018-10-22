//################# TrajXCUDA #########################
//######### Suresh Kondati Natarajan ##################
//##### Lehrstuhl fuer theoretische chemie ############
//######## Ruhr Universitaet Bochum ###################


//life time used by compute_lifetime.cu

#include "cudatools.cuh"

void compute_lifetime(int *A, int *B, float time)
{

  //A is the 0/1 file and B has the number of species to be correlated
   
////////////////////////////////////LIFETIME-PART///////////////////////////////////////////
  
    diffuse_time=time;
  
  
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
	float *Hexch1;

	Hexch1 = (float *) malloc (sizeof(float)*nElt1*SD_store*2);
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
	     }
	 }
	 

int **initHBcol1;// storing the location of the atom in the exch matrix 

initHBcol1 = (int **) malloc (sizeof(int *)*SD_store*2);
     for(int i1=0; i1<SD_store*2; i1++)// for each frame
     {
     initHBcol1[i1] = (int *) malloc (sizeof(int)*B[i1+i*SD_store]);
     for(int j1=0; j1<B[i1+i*SD_store]; j1++)// for number of H bonds this frame
     {
      initHBcol1[i1][j1]=0;
     }
     }
     
   
    for(int i1=0;i1<SD_store*2;i1++)
     {int chker=0;
      for(int j1=0; j1<nElt1; j1++)
       {
	if(A[j1+i1*nElt1+i*SD_store*nElt1] == 1 ) 
        {
	  Hexch1[j1+i1*nElt1]=1.0; //acceptor check
	  initHBcol1[i1][chker] = j1+i1*nElt1; chker+=1;
	}
      }
     }
    
	//cout <<"exch matrix done "<<endl;
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
     //////////////////////////////////////////////////////origin//////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
     for(int j=0; j<SD_store; j=j+skips) //go through restart points or origins
     { 
       if(B[j+i*SD_store] > 0){
       Hexch_sized1 = (float *) malloc (sizeof(float)*B[j+i*SD_store]*SD_store);
       Hexch_sized = (float *) malloc (sizeof(float)*B[j+i*SD_store]*SD_store);
         for(int ia=0;ia<SD_store;ia++)
         {
           for(int ja=0; ja<B[j+i*SD_store]; ja++)////only those present in the first frame in the restart bin
           { 
	        Hexch_sized1[ja+ia*B[j+i*SD_store]] = Hexch1[initHBcol1[j][ja]+ia*nElt1];
	        Hexch_sized[ja+ia*B[j+i*SD_store]] = Hexch1[initHBcol1[j][ja]+ia*nElt1];
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
         for(int ja=0; ja<B[j+i*SD_store]; ja++)///for continuous lifetime
         {
           for(int ia=1;ia<SD_store;ia++)
           {
	     if(Hexch_sized[ja+(ia-1)*B[j+i*SD_store]] == 0) {Hexch_sized[ja+ia*B[j+i*SD_store]] = 0;Hexch_sized1[ja+ia*B[j+i*SD_store]] = 0;} 
	   }
	 }
	 }
	 
	 
	 if(HB_lifestyle=="transient")
	 {
	   //do transient time approximation here
	   
         for(int ja=0; ja<B[j+i*SD_store]; ja++)///for transient intermittent lifetime
         {
           for(int ia=0;ia<SD_store;ia++)
           {
	     if(Hexch_sized[ja+ia*B[j+i*SD_store]] == 0) 
	     {
	       for(int chk=1; chk<=transtime; chk++)
	       {
		 if(ia + chk >= SD_store){break;}
		 else if(Hexch_sized[ja+(ia+chk)*B[j+i*SD_store]] == 1){ Hexch_sized[ja+ia*B[j+i*SD_store]]=1;Hexch_sized1[ja+ia*B[j+i*SD_store]]=1;break;}
	       }
	    } 
	   }
	 }
	 
	 for(int ja=0; ja<B[j+i*SD_store]; ja++)///for imposing continuous lifetime after transient time
         {
           for(int ia=1;ia<SD_store;ia++)
           {
	     if(Hexch_sized[ja+(ia-1)*B[j+i*SD_store]] == 0) {Hexch_sized[ja+ia*B[j+i*SD_store]] = 0;Hexch_sized1[ja+ia*B[j+i*SD_store]] = 0;} 
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
      
      
      
      SD1=(float *) malloc (sizeof(float)*B[j+i*SD_store]*SD_store);
        for(int ia=0;ia<SD_store;ia++)
        {
	 for(int ja=0; ja<B[j+i*SD_store]; ja++)
	 {
           //SD[ja+ia*initHBnum[j+i*SD_store]]=0.0;
	   SD1[ja+ia*B[j+i*SD_store]]=0.0;
	 }
	}
       
      
      dim3 dimBlocka(32,1,32);
      dim3 dimGrida((B[j+i*SD_store]+dimBlocka.x-1)/dimBlocka.x,1,(SD_store+dimBlocka.z-1)/dimBlocka.z);

       
      
      
      float *dev_A1;
            
      cudaMalloc((void **)&dev_A,sizeof(float)*B[j+i*SD_store]*SD_store);
      cudaMalloc((void **)&dev_A1,sizeof(float)*B[j+i*SD_store]*SD_store);
      cudaMalloc((void **)&dev_SD1,sizeof(float)*B[j+i*SD_store]*SD_store);
      
      cudaMemcpy(dev_A,Hexch_sized1,sizeof(float)*B[j+i*SD_store]*SD_store,cudaMemcpyHostToDevice);
      cudaMemcpy(dev_A1,Hexch_sized,sizeof(float)*B[j+i*SD_store]*SD_store,cudaMemcpyHostToDevice);
      cudaMemcpy(dev_SD1,SD1,sizeof(float)*B[j+i*SD_store]*SD_store,cudaMemcpyHostToDevice);

      HBAF_calc<<<dimGrida,dimBlocka>>>(dev_A,dev_A1,dev_SD1,SD_store,B[j+i*SD_store],i,j,origins,skips);

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
      SDreduce<<<dimGrids,dimBlocks>>>(dev_SD1,dev_SDavg1,SD_store,B[j+i*SD_store],fairy);
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
     
 

}
