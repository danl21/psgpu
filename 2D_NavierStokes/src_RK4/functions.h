////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// THIS IS A SOURCE FILE FOR THE 3D_PSPGU CODE                        

// IT CONTAINS HOST FUNCTIONS FOR VARIOUSLY
// WRAPPING FFT ENTRY AND EXIT AND COMPUTING
// STATISTICS AND RECURRENCES FOR FINDING UPO GUESSES
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////

extern "C" void set_gpu_(int *go)
{
  int            dev, deviceCount;
  cudaDeviceProp devProp;

  (cudaGetDeviceCount(&deviceCount));
  if (deviceCount == 0) {
    printf("cutil error: no devices supporting CUDA\n");
    exit(-1);
  }

  (cudaGetDevice(&dev));
  //  (cudaSetDevice());
  (cudaGetDeviceProperties(&devProp,dev));
  printf("\n Using CUDA device %d: %s\n\n", dev,devProp.name);

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////// \
////////////////////////

extern "C" void KR_FFT_ALL(cufftDoubleComplex *d_ZK, cufftDoubleComplex *d_UK,cufftDoubleComplex *d_VK, cufftDoubleReal*d_ZR, cufftDoubleReal*d_UR, cufftDoubleReal*d_VR,cufftHandle PlanZ2D,int *d_ikF){
  cufftDoubleComplex *d_FF;

  //This version does each variable separately (FFT for each Z, U and V rather than one batched FFT)

  (cudaMalloc((void**)&d_FF,sizeof(cufftDoubleComplex)*(ny)*(nx2)));

  conj<<<nblocks,nthreads>>>(d_ZK,d_kx,d_ky);
  setFFN<<<nblocks,nthreads>>>(d_ZK,d_FF,d_ikF);  // Pad truncated wave numbers

  (cufftExecZ2D(PlanZ2D,d_FF,d_ZR));  // Do Z FFT

  conj<<<nblocks,nthreads>>>(d_UK,d_kx,d_ky);
  setFFN<<<nblocks,nthreads>>>(d_UK,d_FF,d_ikF);  // Pad truncated wave numbers

  (cufftExecZ2D(PlanZ2D,d_FF,d_UR));  // Do U FFT

  conj<<<nblocks,nthreads>>>(d_VK,d_kx,d_ky);
  setFFN<<<nblocks,nthreads>>>(d_VK,d_FF,d_ikF);  // Pad truncated wave numbers

  (cufftExecZ2D(PlanZ2D,d_FF,d_VR));  // Do V FFT

  (cudaFree(d_FF));

 }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" double computeNorm(double2 *Zstore,int iStore){
  
  double normZ=0.0;
  for(int i=0; i<nNorm; i++){
    int ik = iNorm[i]+iStore*nkt;
    normZ += Zstore[ik].x*Zstore[ik].x + Zstore[ik].y*Zstore[ik].y;
  }
  return normZ;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" void computeStats(double2 *Z,double2 *ZKtavg,double *kx, double *ky,int *LL){
  
  // function to compute diagnostics; energy, dissipation etc. 
  double eturb = 0.0;
  double dturb = 0.0;
  energy = 0.0;
  diss = 0.0;
  input = 0.0;
  double test = 0.0;

  for(int i=0; i<nkt; i++){
    if(LL[i] ==1 ){
      double wk = max(kx[i]*kx[i] + ky[i]*ky[i],0.001);
      double VZ = Z[i].x*Z[i].x + Z[i].y*Z[i].y;
      diss += VZ;
      energy += VZ/wk;
      ZKtavg[i].x = (ZKtavg[i].x*tavgCount + Z[i].x)/(tavgCount+1.0);
      ZKtavg[i].y = (ZKtavg[i].y*tavgCount + Z[i].y)/(tavgCount+1.0);
      VZ = ZKtavg[i].x*ZKtavg[i].x + ZKtavg[i].y*ZKtavg[i].y;
      eturb += VZ/wk;
      dturb += VZ;
      if(kx[i] == 0.0 && ky[i] == 4.0 ){
        input += - Z[i].x/4.0;
      }
      if(kx[i] == 0.0 && ky[i] == 2.0 ){
        test +=  Z[i].x;
      }
    }
  }
  tavgCount++;
  diss *= 2.0*v2;
  dturb *= 2.0*v2;
  
  fprintf(stats,"%.15e %.15e %.15e %.15e %.15e %.15e %.15e \n",tme,energy,eturb,diss,dturb,input,test);  
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" void computeZdash(double2 *Zstore,
                             double2 *Zdash,
                             double *kx,
                             double *ky,
			     int iStore
                             ){

  // Function to precompute the shifted state vector for recurrence checking
  for(int iShiftY=0; iShiftY < NSHIFTY; iShiftY++){
    for(int iShiftX=0; iShiftX < NSHIFTX; iShiftX++){

      for(int i=0; i<nNorm; i++){
        int ik =iNorm[i];
        int ikk =ik+iStore*nkt;
        int iks = i + nNorm*(iShiftX+NSHIFTX*iShiftY);

        double kkx = kx[ik];
        double kky = ky[ik];
        double shift = kkx*shiftX[iShiftX]+kky*shiftY[iShiftY];

        Zdash[iks].x = cos(shift)*Zstore[ikk].x+sin(shift)*Zstore[ikk].y;
        Zdash[iks].y = cos(shift)*Zstore[ikk].y-sin(shift)*Zstore[ikk].x;
      }
    }
  }
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////

extern "C" void RK_FFT(cufftDoubleComplex *d_ZK, cufftDoubleReal*d_ZR,cufftHandle PlanD2Z,int *d_ikN){
  cufftDoubleComplex *d_FF;

  //Doing a batch of 2 2D FFTs so that now NZK stores the spectral UK and NZK from the old convol in that order.

  (cudaMalloc((void**)&d_FF,sizeof(cufftDoubleComplex)*2*(ny)*(nx2)));

  (cufftExecD2Z(PlanD2Z,d_ZR,d_FF));

  normFF1<<<nblocks,nthreads>>>(d_ZK,d_FF,d_ikN);  // normalise output

  (cudaFree(d_FF));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////

extern "C" void RK_FFT_Z(cufftDoubleComplex *d_ZK, cufftDoubleReal*d_ZR,cufftHandle PlanD2Z,int *d_ikN){
  cufftDoubleComplex *d_FF;

  //Doing a batch of 2 2D FFTs so that now NZK stores the spectral UK and NZK from the old convol in that order.

  (cudaMalloc((void**)&d_FF,sizeof(cufftDoubleComplex)*(ny)*(nx2)));

  (cufftExecD2Z(PlanD2Z,d_ZR,d_FF));

  normFFZ<<<nblocks,nthreads>>>(d_ZK,d_FF,d_ikN);  // normalise output

  (cudaFree(d_FF));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////

extern "C" void find_shift(cufftDoubleComplex *d_Z1, cufftDoubleComplex *d_Z2,int *d_LL,double &Smean){
  double *d_SS,*d_S1,*d_Smean;

  (cudaMalloc((void**)&d_SS,sizeof(double)*(nkt)));
  (cudaMalloc((void**)&d_S1,sizeof(double)*(nkt)));
  (cudaMalloc((void**)&d_Smean,sizeof(double)));

  thrust::device_ptr<double> d_S1ptr(d_S1); // special pointer for global reductions on the GPU
  thrust::device_ptr<double> d_SSptr(d_SS); // special pointer for global reductions on the GPU                        
  // Function to compute the mean shift between two state vectors Z1 and Z2
  mode_shift<<<nblocks,nthreads>>>(d_Z1,d_Z2,d_SS,d_S1,d_kx,d_LL); 
  Smean = thrust::reduce(d_S1ptr,d_S1ptr+nkt,double (0.0), thrust::plus<double>())/double(ikty);

  (cudaMemcpy(d_Smean,&Smean,sizeof(double),cudaMemcpyHostToDevice));
  shift_branch<<<nblocks,nthreads>>>(d_SS,d_S1,d_Smean,d_kx,d_LL);  // normalise output

  Smean = thrust::reduce(d_S1ptr,d_S1ptr+nkt,double (0.0), thrust::plus<double>())/thrust::reduce(d_SSptr,d_SSptr+nkt,double (0.0), thrust::plus<double>());
  
  (cudaFree(d_SS));
  (cudaFree(d_S1));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" void recurrence_check(double2 *Zstore,
				 double *kx,
				 double *ky,
				 int NT,
				 int nCopy
				 ){
  // check for near recurrences
  double minSx,minSy,resPrev,timeUNDER;
  double minRESIDUALinner = 10000.0;
  double minPERIOD = 10000.0;
  double res = 0.0;
  double normFN = 0.0;
  double normFNN = 0.0;
  int minPflag =1;
  int minSTOREinner =0;
  int minCOUNTinner =0;
  double2 Zdash[nNorm*NSHIFTX*NSHIFTY];
  
  computeZdash(Zstore,Zdash,kx,ky,nCopy);
  tme=tstart+(NT-nOut)*delt;
  for(int iCount=(nStore-1); iCount>=1; iCount--){ // note this loop is nStore-1 long
    int iStore;
    double PERIODinner = ((nStore-iCount)*nOut)*(delt);

    // increment over the previous states
    // if nStore = storeSize then we cycle the starting state, iStart
    // this construction below permutes iStore over the storeSize
    iStore = iCount +iStart;
    iStore -= floor(iStore/storeSize)*storeSize;

    //        To save time: when period .gt. 20
    //       only compare evey other previous state.
    if( PERIODinner>20.0 && abs(PERIODinner-freq2*floor(PERIODinner/freq2+0.5))>0.0001)continue;
	  
    resPrev = res;
    res = 10000.0;
    
    for(int iShiftY=0; iShiftY < NSHIFTY; iShiftY++){
      for(int iShiftX=0; iShiftX < NSHIFTX; iShiftX++){
	int is = nNorm*(iShiftX+NSHIFTX*iShiftY);
	double nZdash = 0.0;
	normFNN = 0.0;
	for(int i=0; i<nNorm;i++){
	  int ik = iNorm[i]+iStore*nkt;
	  normFNN -= 2.0*(Zdash[is+i].x*Zstore[ik].x + Zdash[is+i].y*Zstore[ik].y);
	  nZdash += Zdash[is+i].x*Zdash[is+i].x + Zdash[is+i].y*Zdash[is+i].y;
	}
	normFN = normFNN + nZdash + normZStore[iStore];	    
	
	if(sqrt(normFN/normZStore[nCopy]) < res){
	  res = sqrt(normFN/normZStore[nCopy]);
	  minSx = shiftX[iShiftX];
	  minSy = shiftY[iShiftY];
	}
      }
    }
    // Find first turning point as T increases                                                        
    if(res > resPrev && minPflag==1){
      minPERIOD = PERIODinner;
    }else{
      minPflag = 0;
    }
    // If new lower residual is found then store parameters                                           
    if(res < minRESIDUALinner && PERIODinner > minPERIOD){
      minRESIDUALinner = res;
      minSTOREinner =  iStore;
      minCOUNTinner =  nStore - iCount;
      minSHIFTXinner = minSx;
      minSHIFTYinner = minSy;
    }
    
  }
  
  // When min Residual first drops below threshold start new output sequence
  if( minRESIDUALinner < ResidualThreshold && RecOUTflag == 0 ){
    printf("------------------------------------------------\n");
    printf("Time = %e UPO Guess found \n",tme);
    minRESIDUALouter = ResidualThreshold;
    minTIMEinner = tme + minCOUNTinner*nOut*delt;
    RecOUTflag = 1;
  }
  
  // Locate minimum residual for current output sequence and store corresponding data
  if( minRESIDUALinner < minRESIDUALouter && RecOUTflag==1 ){
    minRESIDUALouter = minRESIDUALinner;
    minSTOREouter =  minSTOREinner;
    minTIMEouter = tme + minCOUNTinner*nOut*delt;
    minSHIFTXouter = minSHIFTXinner;
    minSHIFTYouter = minSHIFTYinner;
    RecPERIOD = ((minCOUNTinner)*nOut)*(delt);
  }
  
  // Output best UPO guess at end of current output sequence                                           
  if( minRESIDUALinner > ResidualThreshold && RecOUTflag==1 ){
    RecOUTflag = 0;
    
    // Also compute time spent under residual threshold, normalised by period
    minTIMEinner = tme + minCOUNTinner*nOut*delt- minTIMEinner; 
    timeUNDER = sqrt(minTIMEinner*minTIMEinner/(RecPERIOD*RecPERIOD));

    fwrite(&minTIMEouter,sizeof(double),1,UPOZK);
    fwrite(&v2,sizeof(double),1,UPOZK);
    fwrite(&RecPERIOD,sizeof(double),1,UPOZK);
    fwrite(&minSHIFTXouter,sizeof(double),1,UPOZK);
    fwrite(&minSHIFTYouter,sizeof(double),1,UPOZK);
    fwrite(&minRESIDUALouter,sizeof(double),1,UPOZK);
    fwrite(&timeUNDER,sizeof(double),1,UPOZK);
    fwrite(Zstore+minSTOREouter*nkt,sizeof(double2),nkt,UPOZK);
    izkout = izkout +1;
    printf("Wrote to UPO_Zk.out number %d at time %e \n",izkout,tme);
    
    fprintf(UPOfile,"%d %e %e %e %e %e %e\n",izkout, minTIMEouter,RecPERIOD,minSHIFTXouter,minSHIFTYouter,minRESIDUALouter,timeUNDER);
  }
  fflush(UPOfile);
  fflush(UPOZK);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\
////////////////////////
