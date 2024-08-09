////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// THIS IS A SOURCE FILE FOR THE 3D_PSPGU CODE                                                                                                                                                                                              
// IT CONTAINS HOST FUNCTIONS FOR VARIOUSLY
// WRAPPING FFT ENTRY AND EXIT AND COMPUTING 
// STATISTICS AND RECURRENCES FOR FINDING UPO GUESSES
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
extern "C" void set_gpu_(int *go)
{
  // Function to intialise the GPU very important for the MPI dual
  // instance version
  int            dev, deviceCount;
  cudaDeviceProp devProp;

  (cudaGetDeviceCount(&deviceCount));
  if (deviceCount == 0) {
    printf("cutil error: no devices supporting CUDA\n");
    exit(-1);
  }

  dev = *go;
  printf(" rank = %d \n",dev);
  (cudaSetDevice(*go));
  //  (cudaGetDevice(go));
  (cudaGetDeviceProperties(&devProp,dev));
  printf("\n Using CUDA device %d: %s\n\n", *go,devProp.name);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\
////////////////////////

extern "C" void KR_FFT(cufftDoubleComplex *d_rho, 
		       cufftDoubleReal*d_rx, 
		       cufftHandle PlanZ2D,
		       int *d_ikF,
		       int n2){

  cufftDoubleComplex *d_FF,*d_tmp;

  //This version does each variable separately (FFT for each Z, U and V rather than one batched FFT)

  (cudaMalloc((void**)&d_FF,sizeof(cufftDoubleComplex)*(n2)));
  (cudaMalloc((void**)&d_tmp,sizeof(cufftDoubleComplex)*(nkt)));

  conj<<<nblocks,nthreads>>>(d_rho,d_tmp,d_kx);
  setFFN<<<nblocks,nthreads>>>(d_tmp,d_FF,d_ikF);  // Pad truncated wave numbers
  (cufftExecZ2D(PlanZ2D,d_FF,d_rx));  // Do Z FFT

  (cudaFree(d_FF));
  (cudaFree(d_tmp));
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\
////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\
////////////////////////

extern "C" void KR_FFT_ALL(cufftDoubleComplex *d_xii,
			   cufftDoubleComplex *d_eta,
			   cufftDoubleComplex *d_zet, 
			   cufftDoubleComplex *d_rho, 
			   cufftDoubleComplex *d_UK,
			   cufftDoubleComplex *d_VK,
			   cufftDoubleComplex *d_WK, 
			   cufftDoubleReal*d_xir, 
			   cufftDoubleReal*d_etr, 
			   cufftDoubleReal*d_ztr, 
			   cufftDoubleReal*d_UR, 
			   cufftDoubleReal*d_VR,
			   cufftDoubleReal*d_WR,
			   cufftDoubleReal*d_rx, 
			   cufftDoubleReal*d_ry,
			   cufftDoubleReal*d_rz,
			   cufftHandle PlanZ2D,
			   int *d_ikF, 
			   int n2){

  cufftDoubleComplex *d_FF,*d_tmp;
  // Wrapper for the spectral to physical FFT
  //This version does each variable separately (FFT for each Z, U and V rather than one batched FFT)

  (cudaMalloc((void**)&d_FF,sizeof(cufftDoubleComplex)*(n2)));
  (cudaMalloc((void**)&d_tmp,sizeof(cufftDoubleComplex)*(nkt)));

  conj<<<nblocks,nthreads>>>(d_xii,d_tmp,d_kx);
  setFFN<<<nblocks,nthreads>>>(d_tmp,d_FF,d_ikF);  // Pad truncated wave numbers
  (cufftExecZ2D(PlanZ2D,d_FF,d_xir));  // Do Z FFT

  conj<<<nblocks,nthreads>>>(d_eta,d_tmp,d_kx);
  setFFN<<<nblocks,nthreads>>>(d_tmp,d_FF,d_ikF);  // Pad truncated wave numbers
  (cufftExecZ2D(PlanZ2D,d_FF,d_etr));  // Do Z FFT

  conj<<<nblocks,nthreads>>>(d_zet,d_tmp,d_kx);
  setFFN<<<nblocks,nthreads>>>(d_tmp,d_FF,d_ikF);  // Pad truncated wave numbers
  (cufftExecZ2D(PlanZ2D,d_FF,d_ztr));  // Do Z FFT

  conj<<<nblocks,nthreads>>>(d_VK,d_tmp,d_kx);
  setFFN<<<nblocks,nthreads>>>(d_tmp,d_FF,d_ikF);  // Pad truncated wave numbers
  (cufftExecZ2D(PlanZ2D,d_FF,d_VR));  // Do V FFT

  conj<<<nblocks,nthreads>>>(d_UK,d_tmp,d_kx);
  setFFN<<<nblocks,nthreads>>>(d_tmp,d_FF,d_ikF);  // Pad truncated wave numbers
  (cufftExecZ2D(PlanZ2D,d_FF,d_UR));  // Do U FFT

  conj<<<nblocks,nthreads>>>(d_WK,d_tmp,d_kx);
  setFFN<<<nblocks,nthreads>>>(d_tmp,d_FF,d_ikF);  // Pad truncated wave numbers
  (cufftExecZ2D(PlanZ2D,d_FF,d_WR));  // Do W FFT


  conj<<<nblocks,nthreads>>>(d_rho,d_tmp,d_kx);

  setFFK<<<nblocks,nthreads>>>(d_tmp,d_FF,d_kx,d_ikF);  // Pad truncated wave numbers & do derivative
  (cufftExecZ2D(PlanZ2D,d_FF,d_rx));  // Do rho_x FFT

  setFFK<<<nblocks,nthreads>>>(d_tmp,d_FF,d_ky,d_ikF);  // Pad truncated wave numbers & do derivative
  (cufftExecZ2D(PlanZ2D,d_FF,d_ry));  // Do rho_y FFT

  setFFK<<<nblocks,nthreads>>>(d_tmp,d_FF,d_kz,d_ikF);  // Pad truncated wave numbers & do derivative
  (cufftExecZ2D(PlanZ2D,d_FF,d_rz));  // Do rho_z FFT

  (cudaFree(d_FF));
  (cudaFree(d_tmp));
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\
////////////////////////

extern "C" void RK_FFT(cufftDoubleComplex *d_NXK,
		       cufftDoubleComplex *d_NYK,
		       cufftDoubleComplex *d_NZK, 
		       cufftDoubleComplex *d_RXK, 
		       cufftDoubleReal*d_XR,
		       cufftDoubleReal*d_YR,
		       cufftDoubleReal*d_ZR,
		       cufftDoubleReal*d_rx,		       
		       cufftHandle PlanD2Z,
		       int *d_ikN,
		       int *d_LL,
		       int n2){
  
  cufftDoubleComplex *d_FF;
  // wrapper for physical to spectral FFT
  //Doing a batch of 2 2D FFTs so that now NZK stores the spectral UK and NZK from the old convol in that order.

  (cudaMalloc((void**)&d_FF,sizeof(cufftDoubleComplex)*(n2)));

  (cufftExecD2Z(PlanD2Z,d_XR,d_FF));

  normFF1<<<nblocks,nthreads>>>(d_NXK,d_FF,d_ikN,d_LL);  // normalise output

  (cufftExecD2Z(PlanD2Z,d_YR,d_FF));

  normFF1<<<nblocks,nthreads>>>(d_NYK,d_FF,d_ikN,d_LL);  // normalise output

  (cufftExecD2Z(PlanD2Z,d_ZR,d_FF));

  normFF1<<<nblocks,nthreads>>>(d_NZK,d_FF,d_ikN,d_LL);  // normalise output

  (cufftExecD2Z(PlanD2Z,d_rx,d_FF));

  normFF1<<<nblocks,nthreads>>>(d_RXK,d_FF,d_ikN,d_LL);  // normalise output

  (cudaFree(d_FF));
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\
////////////////////////

extern "C" void RK_FFT_1(cufftDoubleComplex *d_NXK,
			 cufftDoubleReal*d_XR,
			 cufftHandle PlanD2Z,
			 int *d_ikN,
			 int *d_LL,
			 int n2){
  
  cufftDoubleComplex *d_FF;

  // wrapper for single physical to spectral FFT
  (cudaMalloc((void**)&d_FF,sizeof(cufftDoubleComplex)*(n2)));

  (cufftExecD2Z(PlanD2Z,d_XR,d_FF));

  normFF1<<<nblocks,nthreads>>>(d_NXK,d_FF,d_ikN,d_LL);  // normalise output

  (cudaFree(d_FF));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" double computeNorm(double2 *XIIstore,double2 *ETAstore,double2 *ZETstore,double2 *RHOstore,int *LL,int iStore){

  // Compute L2 norm for state vector (used in recurrence check)
  double normX=0.0;
  double normE=0.0;
  double normZ=0.0;
  double normR=0.0;
  for(int i=0; i<nNorm; i++){
    int ik = iNorm[i]+iStore*nkt;
    //    if(LL[iNorm[i]] == 1){ // it may be desirable to use the spectral mask here
	normX += XIIstore[ik].x*XIIstore[ik].x + XIIstore[ik].y*XIIstore[ik].y;
	normE += ETAstore[ik].x*ETAstore[ik].x + ETAstore[ik].y*ETAstore[ik].y;
	normZ += ZETstore[ik].x*ZETstore[ik].x + ZETstore[ik].y*ZETstore[ik].y;
	normR += ri*RHOstore[ik].x*RHOstore[ik].x + RHOstore[ik].y*RHOstore[ik].y;
	//      }
  }
  normZ += normX+normE+normR;
  return normZ;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" void computeBflux(double2 *XIIstore,double2 *ETAstore,double2 *ZETstore,double2 *RHOstore, double *kx, double *ky,double *kz,int *LL,int iStore){

  // !!!!!!!!!!!!!!!!
  // WORK IN PROGRESS
  // !!!!!!!!!!!!!!!!
  // this routine is supposed to classify the scales contributing +ve and -ve buoyancy flux
  double Bhz[nx][nz];
  int cnt[nx][nz];
  int sizeH=0;
  int sizeV=0;
  for(int i=0; i<nx; i++){
    for(int j=0; j<nz; j++){
    Bhz[i][j]=0.0;
    cnt[nx][nz]=0;
    }
  }

  for(int i=0; i<nkt; i++){
    double wkh=1e-10;
    double wkv=1e-10;
    int jh=0;
    int jv=0;
    wkh = max(sqrt(kx[i]*kx[i] + ky[i]*ky[i]),0.001);
    wkv = max(kz[i],0.001);
    jh = int(wkh+0.5);
    jv = int(wkv+0.5);

    int ik = i + iStore*nkt;
    double wk = max(kx[i]*kx[i] + ky[i]*ky[i] + kz[i]*kz[i],0.001);
    double xipR = XIIstore[ik].x;
    double xipC = XIIstore[ik].y;
    double etpR = ETAstore[ik].x;
    double etpC = ETAstore[ik].y;
    double ztpR = ZETstore[ik].x;
    double ztpC = ZETstore[ik].y;
    double ropR = RHOstore[ik].x;
    double ropC = RHOstore[ik].y;

    double vpR = (kz[i]*xipC-kx[i]*ztpC)/wk;
    double vpC = (kx[i]*ztpR-kz[i]*xipR)/wk;
    double wpR = (kx[i]*etpC-ky[i]*xipC)/wk;
    double wpC = (ky[i]*xipR-kx[i]*etpR)/wk;

    double ubpR = vpR*sinTh + wpR*cosTh;
    double ubpC = vpC*sinTh + wpC*cosTh;
    
    double phR = atan2(ropC,ropR);
    double phW = atan2(ubpC,ubpR);

    //    double VZ = ri*sqrt(ubpR*ubpR+ubpC*ubpC)*sqrt(ropR*ropR+ropC*ropC)*cos(phR-phW);
    double VZ = cos(phR-phW);

    Bhz[jh][jv] += VZ;
    cnt[jh][jv]++;
    sizeH = max(sizeH,jh);
    sizeV = max(sizeV,jv);
  }
  for(int i=0; i<sizeH; i++){
    for(int j=0; j<sizeV; j++){
      fprintf(bfile,"%d %d %e\n",i,j,Bhz[i][j]/cnt[i][j]);
    }
    fprintf(bfile,"  \n");
  }
  fflush(bfile);
  
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" void computeStats(double2 *XIIstore,double2 *ETAstore,double2 *ZETstore,double2 *RHOstore,double2 *XIItavg,double2 *ETAtavg,double2 *ZETtavg, double2 *RHOtavg, double *kx, double *ky,double *kz,int *LL,int iStore){

  // Main routine for computing diagnostic quantities, particularly energetic ones
  double eturb = 0.0;
  double dturb = 0.0;
  double peturb = 0.0;
  double pedturb = 0.0;
  double ped1 =0.0;
  double ped2 =0.0;

  energy = 0.0;
  diss = 0.0;
  input = 0.0;
  bflux = 0.0;
  pe = 0.0;
  ped = 0.0;
  
  for(int i=0; i<nkt; i++){
    if(LL[i] == 1){
    int ik = i + iStore*nkt;
    double wk = max(kx[i]*kx[i] + ky[i]*ky[i]+kz[i]*kz[i],0.001);

    double VZX = XIIstore[ik].x*XIIstore[ik].x + XIIstore[ik].y*XIIstore[ik].y;
    double VZY = ETAstore[ik].x*ETAstore[ik].x + ETAstore[ik].y*ETAstore[ik].y;
    double VZZ = ZETstore[ik].x*ZETstore[ik].x + ZETstore[ik].y*ZETstore[ik].y;
    double VZ = VZX+VZY+VZZ;
    diss += VZ; // total kinetic energy dissipation <|grad u|^2>
    energy += VZ/wk; // total kinetic energy <u^2>

    if(kx[i] ==0.0 && ky[i] == kfy && kz[i] == 0.0){
      input -= ampfor*ZETstore[ik].x/kfy; // energy input <u.f>
    }

    double VR = RHOstore[ik].x*RHOstore[ik].x + RHOstore[ik].y*RHOstore[ik].y;
    pe += VR; // total density variance <rho^2>
    ped+= VR*wk; // total density variance dissipation (chi) <|grad(rho)|^2>
    ped1+= kx[i]*kx[i]*VR; // chi contributions from various gradients
    ped2+= ky[i]*ky[i]*VR;

    // compute running mean for the purposes of a running estimate of the fluctuations
    XIItavg[i].x = (XIItavg[i].x*tavgCount + XIIstore[ik].x)/(tavgCount+1.0);
    XIItavg[i].y = (XIItavg[i].y*tavgCount + XIIstore[ik].y)/(tavgCount+1.0);
    ETAtavg[i].x = (ETAtavg[i].x*tavgCount + ETAstore[ik].x)/(tavgCount+1.0);
    ETAtavg[i].y = (ETAtavg[i].y*tavgCount + ETAstore[ik].y)/(tavgCount+1.0);
    ZETtavg[i].x = (ZETtavg[i].x*tavgCount + ZETstore[ik].x)/(tavgCount+1.0);
    ZETtavg[i].y = (ZETtavg[i].y*tavgCount + ZETstore[ik].y)/(tavgCount+1.0);
    RHOtavg[i].x = (RHOtavg[i].x*tavgCount + RHOstore[ik].x)/(tavgCount+1.0);
    RHOtavg[i].y = (RHOtavg[i].y*tavgCount + RHOstore[ik].y)/(tavgCount+1.0);

    VZX = XIItavg[i].x*XIItavg[i].x + XIItavg[i].y*XIItavg[i].y;
    VZY = ETAtavg[i].x*ETAtavg[i].x + ETAtavg[i].y*ETAtavg[i].y;
    VZZ = ZETtavg[i].x*ZETtavg[i].x + ZETtavg[i].y*ZETtavg[i].y;
    VZ = VZX+VZY+VZZ;

    VR = RHOtavg[i].x*RHOtavg[i].x + RHOtavg[i].y*RHOtavg[i].y;

    // compute turbulent energetic quantities
    eturb += VZ/wk;
    dturb += VZ;
    peturb += VR;
    pedturb += VR*wk;

    // put various components in local variables (efficient for cache use)
    double xipR = XIIstore[ik].x;
    double xipC = XIIstore[ik].y;
    double etpR = ETAstore[ik].x;
    double etpC = ETAstore[ik].y;
    double ztpR = ZETstore[ik].x;
    double ztpC = ZETstore[ik].y;
    double ropR = RHOstore[ik].x;
    double ropC = RHOstore[ik].y;
 
    // velocity components in spectral space
    double vpR = (kz[i]*xipC-kx[i]*ztpC)/wk;
    double vpC = (kx[i]*ztpR-kz[i]*xipR)/wk;
    double wpR = (kx[i]*etpC-ky[i]*xipC)/wk;
    double wpC = (ky[i]*xipR-kx[i]*etpR)/wk;

    double ubpR = vpR*sinTh + wpR*cosTh;
    double ubpC = vpC*sinTh + wpC*cosTh;
    // phase between product terms in buoyancy flux
    double phR = atan2(ropC,ropR);
    double phW = atan2(ubpC,ubpR);
    // buoyancy flux <rho*q>
    bflux += ri*sqrt(ubpR*ubpR+ubpC*ubpC)*sqrt(ropR*ropR+ropC*ropC)*cos(phR-phW);

    }
  }
  // *****************
  // Here we update the forcing amplitude using the throttling method
  // designed to maintain a target dissipation rate and stationary energetics
  if(dtarget !=0.0){
    ampfor = ampfor*(dtarget-bflux)/input;
  }
  // *****************

  // mean quantities
  dissAvg = (dissAvg*tavgCount + diss)/(tavgCount+1.0);
  enerAvg = (enerAvg*tavgCount + energy)/(tavgCount+1.0);
  peAvg = (peAvg*tavgCount + pe)/(tavgCount+1.0);
  pedAvg = (pedAvg*tavgCount + ped)/(tavgCount+1.0);
  inputAvg = (inputAvg*tavgCount + input)/(tavgCount+1.0);
  bfluxAvg = (bfluxAvg*tavgCount + bflux)/(tavgCount+1.0);
  mixAvg = ((ped/(diss/ri+ped))*tavgCount + mixAvg)/(tavgCount+1.0);

  // multiply by various parameters (N.B. output D/Dlam and I/Dlam)
  input *= v2*kfy*kfy*2.0;
  diss *= v2*v2*4.0*kfy*kfy;
  dturb = dissAvg-dturb;
  dturb *= v2*v2*4.0*kfy*kfy;
  eturb = enerAvg-eturb;
  peturb = peAvg-peturb;
  pedturb= pedAvg-pedturb;
    
  dturbAvg = (dturbAvg*tavgCount + dturb)/(tavgCount+1.0);
  eturbAvg = (eturbAvg*tavgCount + eturb)/(tavgCount+1.0);
  //now dissipation stats normalised with laminar

  dissMax = max(diss,dissMax);
  dissMin = min(diss,dissMin);
  enerMax = max(energy,enerMax);
  enerMin = min(energy,enerMin);
  inputMax = max(input,inputMax);
  inputMin = min(input,inputMin);

  // store history for computing p.d.f. at the end
  Diss[countStats-1]  = diss;
  Energy[countStats-1]= energy;
  Input[countStats-1] = input;

  tavgCount++;
  fprintf(stats,"%e %e %e %e %e %e %e \n",tme,energy,eturb,diss,dturb,input,Umax);
  fflush(stats);
  fprintf(engy,"%e %e %e %e %e %e %e %e\n",tme,bflux,pe,peturb,ped,pedturb,ped1,ped2);
  fflush(engy);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

extern "C" void computeZdash(double2 *XIIstore,
			     double2 *ETAstore,
			     double2 *ZETstore,
			     double2 *RHOstore,
			     double2 *XIIdash,
			     double2 *ETAdash,
                             double2 *ZETdash,
                             double2 *RHOdash,
                             double *kx,
                             double *ky,
                             double *kz,
                             int iStore
                             ){

  // Store the shifted components of the current state vector
  // This is to reduce the workload on every check with the
  // history during recurrence search. 

  // first loop over the shifts
  for(int iShiftY=0; iShiftY < NSHIFTY; iShiftY++){
    for(int iShiftZ=0; iShiftZ < NSHIFTZ; iShiftZ++){
      for(int iShiftX=0; iShiftX < NSHIFTX; iShiftX++){
      
	for(int i=0; i<nNorm; i++){
	  int ik =iNorm[i];
	  int ikk =ik+iStore*nkt;
	  int iks = i + nNorm*(iShiftX+NSHIFTX*(iShiftY+NSHIFTY*iShiftZ));

	  double kkx = kx[ik];
	  double kky = ky[ik];
	  double kkz = kz[ik];
	  double shift = kkx*shiftX[iShiftX]+kky*shiftY[iShiftY]+kkz*shiftZ[iShiftZ]; // compute the shift multiplier

	  // shift the components and store in the *dash variables
	  XIIdash[iks].x = cos(shift)*XIIstore[ikk].x+sin(shift)*XIIstore[ikk].y;
	  XIIdash[iks].y = cos(shift)*XIIstore[ikk].y-sin(shift)*XIIstore[ikk].x;

	  ETAdash[iks].x = cos(shift)*ETAstore[ikk].x+sin(shift)*ETAstore[ikk].y;
	  ETAdash[iks].y = cos(shift)*ETAstore[ikk].y-sin(shift)*ETAstore[ikk].x;

	  ZETdash[iks].x = cos(shift)*ZETstore[ikk].x+sin(shift)*ZETstore[ikk].y;
	  ZETdash[iks].y = cos(shift)*ZETstore[ikk].y-sin(shift)*ZETstore[ikk].x;

	  RHOdash[iks].x = cos(shift)*RHOstore[ikk].x+sin(shift)*RHOstore[ikk].y;
	  RHOdash[iks].y = cos(shift)*RHOstore[ikk].y-sin(shift)*RHOstore[ikk].x;
	}
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" void recurrence_check(double2 *XIIstore,
				 double2 *ETAstore,
				 double2 *ZETstore,
				 double2 *RHOstore,
                                 double *kx,
                                 double *ky,
                                 double *kz,
                                 int NT,
                                 int nCopy
                                 ){

  // This function loops over the history in *store variables and compares 
  // agains the shifted current stated vector for near recurrence. 
  // There then follows a bit of complicated logic to find the smallest 
  // residual for a given recurrent episode. (I think this can be improved- for future work)

  double minSx,minSy,minSz,resPrev;
  double minRESIDUALinner = 10000.0;
  double minPERIOD = 10000.0;
  double res = 0.0;
  double normFN = 0.0;
  double normFNN = 0.0;
  int minPflag =1;
  int minSTOREinner =0;
  int minCOUNTinner =0;
  double2 *XIIdash,*ETAdash,*ZETdash,*RHOdash;

  XIIdash = (double2*)malloc(sizeof(double2)*nNorm*NSHIFTX*NSHIFTY*NSHIFTZ);
  ETAdash = (double2*)malloc(sizeof(double2)*nNorm*NSHIFTX*NSHIFTY*NSHIFTZ);
  ZETdash = (double2*)malloc(sizeof(double2)*nNorm*NSHIFTX*NSHIFTY*NSHIFTZ);
  RHOdash = (double2*)malloc(sizeof(double2)*nNorm*NSHIFTX*NSHIFTY*NSHIFTZ);
  
  // get the shifted current state (see above)
  computeZdash(XIIstore,ETAstore,ZETstore,RHOstore,XIIdash,ETAdash,ZETdash,RHOdash,kx,ky,kz,nCopy);

  tme=tstart+(NT-nOut)*delt;
  // start the loop through the history, starting at the nearest and going backwards
  for(int iCount=(nStore-1); iCount>=1; iCount--){ // note this loop is nStore-1 long
    int iStore;
    double PERIODinner = ((nStore-iCount)*nOut)*(delt); // period given by the historical proximity

    // increment over the previous states
    // if nStore = storeSize then we cycle the starting state, iStart
    // this construction below permutes iStore over the storeSize
    iStore = iCount +iStart;
    iStore -= floor(iStore/storeSize)*storeSize;

    //        To save time: when period .gt. 20
    //       only compare evey other previous state.
    if(PERIODinner<0.5)continue;
    if( PERIODinner>20.0 && abs(PERIODinner-freq2*floor(PERIODinner/freq2+0.5))>0.0001)continue;

    resPrev = res;
    res = 10000.0;

    for(int iShiftY=0; iShiftY < NSHIFTY; iShiftY++){
      for(int iShiftZ=0; iShiftZ < NSHIFTZ; iShiftZ++){
	for(int iShiftX=0; iShiftX < NSHIFTX; iShiftX++){
	  int is = nNorm*(iShiftX+NSHIFTX*(iShiftY+NSHIFTY*iShiftZ));
	  double nZdash = 0.0;
	  normFNN = 0.0;
	  // compute the L2 form of the difference between the shifted current vector and the historical one
	  for(int i=0; i<nNorm;i++){
	    int ik = iNorm[i]+iStore*nkt;

	    normFNN -= 2.0*(XIIdash[is+i].x*XIIstore[ik].x + XIIdash[is+i].y*XIIstore[ik].y);
	    nZdash += XIIdash[is+i].x*XIIdash[is+i].x + XIIdash[is+i].y*XIIdash[is+i].y;

	    normFNN -= 2.0*(ETAdash[is+i].x*ETAstore[ik].x + ETAdash[is+i].y*ETAstore[ik].y);
	    nZdash += ETAdash[is+i].x*ETAdash[is+i].x + ETAdash[is+i].y*ETAdash[is+i].y;

	    normFNN -= 2.0*(ZETdash[is+i].x*ZETstore[ik].x + ZETdash[is+i].y*ZETstore[ik].y);
	    nZdash += ZETdash[is+i].x*ZETdash[is+i].x + ZETdash[is+i].y*ZETdash[is+i].y;

	    normFNN -= ri*(2.0*(RHOdash[is+i].x*RHOstore[ik].x + RHOdash[is+i].y*RHOstore[ik].y));
	    nZdash += ri*(RHOdash[is+i].x*RHOdash[is+i].x + RHOdash[is+i].y*RHOdash[is+i].y);
	  }
        normFN = normFNN + nZdash + normZStore[iStore];

	// find the smallest residual across this sweep and store the shifts
	if(sqrt(normFN/(normZStore[iStore]+normZStore[nCopy])) < res){
	  res = sqrt(normFN/(normZStore[iStore]+normZStore[nCopy]));
          minSx = shiftX[iShiftX];
          minSy = shiftY[iShiftY];
          minSz = shiftZ[iShiftZ];
        }
	}
      }
    }
    // Find first turning point of the residual as T increases
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
      minSHIFTZinner = minSz;
    }

  }
	
  // When min Residual first drops below threshold start new output sequence
  if( minRESIDUALinner < ResidualThreshold && RecOUTflag == 0 ){ 
    printf("------------------------------------------------\n");
    printf("Time = %e UPO Guess found \n",tme);
    minRESIDUALouter = ResidualThreshold;
    RecOUTflag = 1;
  }

  // Locate minimum residual for current output sequence and store corresponding data
  if( minRESIDUALinner < minRESIDUALouter && RecOUTflag==1 ){
    minRESIDUALouter = minRESIDUALinner;
    minSTOREouter =  minSTOREinner;
    minTIMEouter = tme + minCOUNTinner*nOut*delt;
    minSHIFTXouter = minSHIFTXinner;
    minSHIFTYouter = minSHIFTYinner;
    minSHIFTZouter = minSHIFTZinner;
    RecPERIOD = ((minCOUNTinner)*nOut)*(delt);
  }

  // Output best UPO guess at end of current output sequence
  if( minRESIDUALinner > ResidualThreshold && RecOUTflag==1 ){
    RecOUTflag = 0;

    fwrite(&minTIMEouter,sizeof(double),1,UPOZK);
    fwrite(&v2,sizeof(double),1,UPOZK);
    fwrite(&RecPERIOD,sizeof(double),1,UPOZK);
    fwrite(&minSHIFTXouter,sizeof(double),1,UPOZK);
    fwrite(&minSHIFTYouter,sizeof(double),1,UPOZK);
    fwrite(&minSHIFTZouter,sizeof(double),1,UPOZK);
    fwrite(&minRESIDUALouter,sizeof(double),1,UPOZK);
    fwrite(XIIstore+minSTOREouter*nkt,sizeof(double2),nkt,UPOZK);
    fwrite(ETAstore+minSTOREouter*nkt,sizeof(double2),nkt,UPOZK);
    fwrite(ZETstore+minSTOREouter*nkt,sizeof(double2),nkt,UPOZK);
    fwrite(RHOstore+minSTOREouter*nkt,sizeof(double2),nkt,UPOZK);
    izkout = izkout +1;
    printf("Wrote to UPO_Zk.out number %d at time %e \n",izkout,tme);

    fprintf(UPOfile,"%d %e %e %e %e %e %e \n",izkout, minTIMEouter,RecPERIOD,minSHIFTXouter,minSHIFTYouter,minSHIFTZouter,minRESIDUALouter);
  }

  free(XIIdash);
  free(ETAdash);
  free(ZETdash);
  free(RHOdash);
  fflush(UPOfile);
  fflush(UPOZK);
}


