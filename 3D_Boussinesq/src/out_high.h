// ------------------------------------------------------------------
// high frequency outputs and recurrence checking
// ------------------------------------------------------------------

    if(tcheck < 0.0 && (*RecFLAG ==1 || *statsFLAG!=0) && timestep!=0){
      countStats++;

      if(nStore >= 1 && *RecFLAG==1){
        recurrence_check(XIIstore,ETAstore,ZETstore,RHOstore,kx,ky,kz,timestep,nCopy); // check for near recurrence
        fflush(stdout);
      }

      // keep track of indexing for the copy index and the back search start index
      // copy index should lag the start index
      nStore++;
      if(nStore >= storeSize){
        nStore = storeSize;
        nCopy++;
        nCopy -= floor(nCopy/storeSize)*storeSize;
        iStart=nCopy;
      }else{
        nCopy = nStore;
      }
      cudaThreadSynchronize();

      // copy historical record from the GPU
      (cudaMemcpy(XIIstore+nCopy*nkt,d_xii,sizeof(cufftDoubleComplex)*nkt,cudaMemcpyDeviceToHost));
      (cudaMemcpy(ETAstore+nCopy*nkt,d_eta,sizeof(cufftDoubleComplex)*nkt,cudaMemcpyDeviceToHost));
      (cudaMemcpy(ZETstore+nCopy*nkt,d_zet,sizeof(cufftDoubleComplex)*nkt,cudaMemcpyDeviceToHost));
      (cudaMemcpy(RHOstore+nCopy*nkt,d_rho,sizeof(cufftDoubleComplex)*nkt,cudaMemcpyDeviceToHost));

      if(*RecFLAG==1)normZStore[nCopy]=computeNorm(XIIstore,ETAstore,ZETstore,RHOstore,LC,nCopy); // compute the current ||Z||

      if(*statsFLAG==1)computeStats(XIIstore,ETAstore,ZETstore,RHOstore,XIItavg,ETAtavg,ZETtavg,RHOtavg,kx,ky,kz,LC,nCopy);  // get the stats    

      if(*statsFLAG==2 && tme >50.0 )computeStats(XIIstore,ETAstore,ZETstore,RHOstore,XIItavg,ETAtavg,ZETtavg,RHOtavg,kx,ky,kz,LC,nCopy); // get the stats at t>T_1

      // We might have not allocated enough space for the diagnostics, reallocate it... (needed for pdfs)
      if(countStats>nStats && *statsFLAG !=0){
	nStats=2*nStats;
	double *tempDiss=(double*)realloc(Diss,sizeof(double)*nStats);
	double *tempEnergy=(double*)realloc(Energy,sizeof(double)*nStats);
	double *tempInput=(double*)realloc(Input,sizeof(double)*nStats);
	if( tempDiss != NULL){
	  Diss = tempDiss;
	}else{
	  printf("Error reallocating Diss, stats before here... %e ",tme);
	  countStats=0;
	}
	if( tempEnergy != NULL){
	  Energy = tempEnergy;
	}else{
	  printf("Error reallocating Energy, stats before here... %e ",tme);
	  countStats=0;
	}
	if( tempInput != NULL){
	  Input = tempInput;
	}else{
	  printf("Error reallocating Input, stats before here... %e ",tme);
	  countStats=0;
	}
      }
    }
// ------------------------------------------------------------------
// ------------------------------------------------------------------
// ------------------------------------------------------------------
