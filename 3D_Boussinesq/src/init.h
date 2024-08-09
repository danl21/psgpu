int ikty=*IKTY;

// Define CPU variables                                                                                                                                                                                                                    
ny = *NY;
nx = *NX;
nz = *NZ;
kfy= *KFY;
nkt = *NKT;
nr  = (*NY)*(*NX)*(*NZ);
nx2 = (*NX)/2 +1;
n2 = nx2*(*NY)*(*NZ);
in = 1.0/nr;
dtarget = *Dtarget;
ri = *Ri;
nOut = *NOUT;
tme = *TSTART;
tstart = *TSTART;
ampfor = -*AMPFOR;
v2 = *V2;
ResidualThreshold = *ResThresh;
rank = *RANK;

// set number of timesteps/timestep                                                                                                                                                                                                        
double nStop = *Tstep/(*DELT);
if(nStop**DELT != *Tstep){
  delt = *Tstep/nStop;
 }else{
  delt = *DELT;
 }

double C = twopi/nx; //Courant number                                                                                                                                                                                                      
double Tout = nOut*delt; // total T                                                                                                                                                                                                        
double Vout = (*NOUTV)*delt; // output interval                                                                                                                                                                                            
double deltp=delt; // max timestep when adaptivity is on                                                                                                                                                                                   

// write params to stdout for checking
printf("Tstep %e \n",*Tstep);
printf("ampfor %e \n",ampfor);
printf("delt %e \n",delt);
printf("v2 %e \n",v2);
printf("Ri %e\n",*Ri);
printf("Sc %e\n",*Sc);
printf("theta %e \n",*Theta);
printf("kfy %e \n",kfy);

  // intialising various bulk measures for averaging or min/max
  double IN = 1.0/nr;
  double RiGMax = 0.0;
  dissMax = 0.0;
  dissMin = 10000.0;
  inputMax = 0.0;
  inputMin = 10000.0;
  enerMax = 0.0;
  enerMin = 10000.0;
  enerAvg = 0.0;
  dissAvg = 0.0;
  pedAvg = 0.0;
  mixAvg  = 0.0;
  inputAvg = 0.0;
  eturbAvg = 0.0;
  dturbAvg = 0.0;
  bfluxAvg = 0.0;
  nStore = 0;
  freq2 = 2*nOut*delt;
  tavgCount =0;
  sinTh = sin(*Theta); // if gravity is inclined we need the sine and cosine of the angle for projecting
  cosTh = cos(*Theta);

  // Various host arrays for storing history and averaging/writing to disc.
  int nCopy;
  double2 *XIIstore,*ETAstore,*ZETstore,*RHOstore;
  double2 *XIItavg,*ETAtavg,*ZETtavg,*RHOtavg;
  double  *Umean,*Vmean,*Wmean,*uhat,*vhat,*what,*ztr,*xir,*etr,*ro;
  size_t avail,total;
  
  // Define global device variables
  int *d_ikF, *d_ikN, *d_LL; // arrays for indexing and masking
  cufftDoubleComplex *d_xii,*d_eta,*d_zet,*d_rho,*d_UK, *d_VK, *d_WK,*d_rk; // basic Fourier components
  cufftDoubleComplex *d_x0,*d_x1,*d_x2,*d_x3,*d_e0,*d_e1,*d_e2,*d_e3,*d_z0,*d_z1,*d_z2,*d_z3,*d_r0,*d_r1,*d_r2,*d_r3; // RK4 arrays
  cufftDoubleReal *d_xir,*d_etr,*d_ztr,*d_rx,*d_ry,*d_rz,*d_UR,*d_VR,*d_WR,*d_ro; // basic physical components
  double *d_Urms,*d_Vrms,*d_Wrms; // rms velocities
  cufftHandle PlanZ2D,PlanD2Z; // FFT plans

  // Allocate global memory on GPU. (Constant memory does not need allocating) 	
  (cudaMalloc((void**)&d_xii,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_x0,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_x1,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_x2,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_x3,sizeof(cufftDoubleComplex)*(nkt)));

  (cudaMalloc((void**)&d_eta,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_e0,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_e1,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_e2,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_e3,sizeof(cufftDoubleComplex)*(nkt)));

  (cudaMalloc((void**)&d_zet,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_z0,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_z1,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_z2,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_z3,sizeof(cufftDoubleComplex)*(nkt)));


  (cudaMalloc((void**)&d_rho,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_r0,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_r1,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_r2,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_r3,sizeof(cufftDoubleComplex)*(nkt)));

  (cudaMalloc((void**)&d_UK,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_VK,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_WK,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_rk,sizeof(cufftDoubleComplex)*(nkt)));
  
  (cudaMalloc((void**)&d_xir,sizeof(cufftDoubleReal)*(nr)));
  (cudaMalloc((void**)&d_etr,sizeof(cufftDoubleReal)*(nr)));
  (cudaMalloc((void**)&d_ztr,sizeof(cufftDoubleReal)*(nr)));
  (cudaMalloc((void**)&d_rx,sizeof(cufftDoubleReal)*(nr)));
  (cudaMalloc((void**)&d_ry,sizeof(cufftDoubleReal)*(nr)));
  (cudaMalloc((void**)&d_rz,sizeof(cufftDoubleReal)*(nr)));
  (cudaMalloc((void**)&d_ro,sizeof(cufftDoubleReal)*(nr)));
  (cudaMalloc((void**)&d_UR,sizeof(cufftDoubleReal)*nr));
  (cudaMalloc((void**)&d_VR,sizeof(cufftDoubleReal)*nr));
  (cudaMalloc((void**)&d_WR,sizeof(cufftDoubleReal)*nr));
  
  (cudaMalloc((void**)&d_kx,sizeof(double)*(nkt)));
  (cudaMalloc((void**)&d_ky,sizeof(double)*(nkt)));
  (cudaMalloc((void**)&d_kz,sizeof(double)*(nkt)));
  (cudaMalloc((void**)&d_LL,sizeof(int)*nkt));
  (cudaMalloc((void**)&d_ikN,sizeof(int)*nkt));
  (cudaMalloc((void**)&d_ikF,sizeof(int)*n2));
  
  // Copy state data to GPU global memory 
  (cudaMemcpy(d_xii,xii,sizeof(cufftDoubleComplex)*(nkt),cudaMemcpyHostToDevice));
  (cudaMemcpy(d_eta,eta,sizeof(cufftDoubleComplex)*(nkt),cudaMemcpyHostToDevice));
  (cudaMemcpy(d_zet,zet,sizeof(cufftDoubleComplex)*(nkt),cudaMemcpyHostToDevice));
  (cudaMemcpy(d_rho,rho,sizeof(cufftDoubleComplex)*(nkt),cudaMemcpyHostToDevice));
  //Set up various arrays to enable generic kernel calls
  //i.e. calculate indexing for padding either side of FFTs, wavenumber arrays, mask, and timestep arrays.
  // This must be done on CPU for scalability (large problems violate max threads per block)
  (cudaMemcpy(d_kx,kx,sizeof(double)*(nkt),cudaMemcpyHostToDevice));
  (cudaMemcpy(d_ky,ky,sizeof(double)*(nkt),cudaMemcpyHostToDevice));
  (cudaMemcpy(d_kz,kz,sizeof(double)*(nkt),cudaMemcpyHostToDevice));
  (cudaMemcpy(d_ikF,ikF, sizeof(int)*n2,cudaMemcpyHostToDevice));
  (cudaMemcpy(d_ikN,ikN,sizeof(int)*nkt,cudaMemcpyHostToDevice));
  (cudaMemcpy(d_LL,LC,sizeof(int)*nkt,cudaMemcpyHostToDevice));

  // Copy constsant parameters to GPU constant memory
  (cudaMemcpyToSymbol(d_IN,&IN,sizeof(double)));
  (cudaMemcpyToSymbol(d_AMPFOR,AMPFOR,sizeof(double)));
  (cudaMemcpyToSymbol(d_DELT,&delt,sizeof(double)));
  (cudaMemcpyToSymbol(d_v2,V2,sizeof(double)));
  (cudaMemcpyToSymbol(d_Ri,Ri,sizeof(double)));
  (cudaMemcpyToSymbol(d_Sc,Sc,sizeof(double))); 
  (cudaMemcpyToSymbol(d_SinTh,&sinTh,sizeof(double)));
  (cudaMemcpyToSymbol(d_CosTh,&cosTh,sizeof(double)));
  (cudaMemcpyToSymbol(d_KFY,KFY,sizeof(double)));
  (cudaMemcpyToSymbol(d_IKTX,IKTX,sizeof(int)));
  (cudaMemcpyToSymbol(d_IKTY,IKTY,sizeof(int)));
  (cudaMemcpyToSymbol(d_IKTZ,IKTZ,sizeof(int)));
  (cudaMemcpyToSymbol(d_KTZ,KTZ,sizeof(int)));
  (cudaMemcpyToSymbol(d_NX,NX,sizeof(int)));
  (cudaMemcpyToSymbol(d_NX2,&nx2,sizeof(int)));
  (cudaMemcpyToSymbol(d_NY,NY,sizeof(int)));
  (cudaMemcpyToSymbol(d_NZ,NY,sizeof(int)));
  (cudaMemcpyToSymbol(d_OR,&nr,sizeof(int)));
  (cudaMemcpyToSymbol(d_OK,&nkt,sizeof(int)));
  (cudaMemcpyToSymbol(d_O2,&n2,sizeof(int)));
    
  //Set FFT Plans
  (cufftPlan3d(&PlanZ2D,*NZ,*NY,*NX,CUFFT_Z2D));
  (cufftPlan3d(&PlanD2Z,*NZ,*NY,*NX,CUFFT_D2Z));
  thrust::device_ptr<double> d_Uptr(d_UR); // special variable for global reductions on the GPU
  //allocate real arrays for writing
  xir=(double*)malloc(sizeof(double)*nr);
  ztr=(double*)malloc(sizeof(double)*nr);
  etr=(double*)malloc(sizeof(double)*nr);
  ro =(double*)malloc(sizeof(double)*nr);

  // set the output file names
  char xiifile[sizeof "JOB#/xii000.dat"];
  char etafile[sizeof "JOB#/eta000.dat"];
  char zetfile[sizeof "JOB#/zet000.dat"];
  char rhofile[sizeof "JOB#/rho000.dat"];
  int nStats = 100**Tstep/Tout+1;
  printf("nStats= %d \n",nStats);
  FILE* sweep;
  // open the output streams
  if(*statsFLAG !=0){    
    sweep = fopen("sweep.dat","a");

    char statsfile[sizeof "JOB#/stats.dat"];
    sprintf(statsfile, "JOB%1d/stats.dat",rank);
    stats = fopen(statsfile,"w");

    char engyfile[sizeof "JOB#/engy.dat"];
    sprintf(engyfile, "JOB%1d/engy.dat",rank);
    engy = fopen(engyfile,"w");

    char avgsfile[sizeof "JOB#/avgs.dat"];
    sprintf(avgsfile, "JOB%1d/avgs.dat",rank);
    avgs = fopen(avgsfile,"w");

    char meansfile[sizeof "JOB#/means.dat"];
    sprintf(meansfile, "JOB%1d/means.dat",rank);
    means = fopen(meansfile,"w");

    char pdfsfile[sizeof "JOB#/pdfs.dat"];
    sprintf(pdfsfile, "JOB%1d/pdfs.dat",rank);
    pdfs = fopen(pdfsfile,"w");

    char specfile[sizeof "JOB#/spec.dat"];
    sprintf(specfile, "JOB%1d/spec.dat",rank);
    spec = fopen(specfile,"w");

    char RiGfile[sizeof "JOB#/RiG.dat"];
    sprintf(RiGfile, "JOB%1d/RiG.dat",rank);
    RiG = fopen(RiGfile,"w");

    // Allocate the host arrays for history and averages
    XIIstore=(double2*)malloc(sizeof(double2)*nkt*storeSize);
    ETAstore=(double2*)malloc(sizeof(double2)*nkt*storeSize);
    ZETstore=(double2*)malloc(sizeof(double2)*nkt*storeSize);
    RHOstore=(double2*)malloc(sizeof(double2)*nkt*storeSize);
    XIItavg=(double2*)malloc(sizeof(double2)*nkt);
    ETAtavg=(double2*)malloc(sizeof(double2)*nkt);
    ZETtavg=(double2*)malloc(sizeof(double2)*nkt);
    RHOtavg=(double2*)malloc(sizeof(double2)*nkt);

    // Allocate the device arrays for rmses
    (cudaMalloc((void**)&d_Urms,sizeof(cufftDoubleReal)*nr));
    (cudaMalloc((void**)&d_Vrms,sizeof(cufftDoubleReal)*nr));
    (cudaMalloc((void**)&d_Wrms,sizeof(cufftDoubleReal)*nr));

    // Allocate host arrays for energetic history
    Diss=(double*)malloc(sizeof(double)*nStats);
    Energy=(double*)malloc(sizeof(double)*nStats);
    Input=(double*)malloc(sizeof(double)*nStats);

    // initalise to zero (perhaps unecessary?)
    zeroReal<<<nblocks,nthreads>>>(d_Urms);
    zeroReal<<<nblocks,nthreads>>>(d_Vrms);
    zeroReal<<<nblocks,nthreads>>>(d_Wrms);

    for(int i=0; i<nkt; i++){
      XIItavg[i].x=0.0;
      XIItavg[i].y=0.0;

      ETAtavg[i].x=0.0;
      ETAtavg[i].y=0.0;

      ZETtavg[i].x=0.0;
      ZETtavg[i].y=0.0;
    }

    iStart = 0;
    nCopy = 0;
    // Copy I.C. to history
    (cudaMemcpy(XIIstore,d_xii,sizeof(cufftDoubleComplex)*nkt,cudaMemcpyDeviceToHost));
    (cudaMemcpy(ETAstore,d_eta,sizeof(cufftDoubleComplex)*nkt,cudaMemcpyDeviceToHost));
    (cudaMemcpy(ZETstore,d_zet,sizeof(cufftDoubleComplex)*nkt,cudaMemcpyDeviceToHost));
    (cudaMemcpy(RHOstore,d_rho,sizeof(cufftDoubleComplex)*nkt,cudaMemcpyDeviceToHost));
    
  }
  // Do a check of global memory use
  avail =0;
  total = 0;
  cudaMemGetInfo(&avail,&total);
  
  printf("\n total : %f MB \n",float(total)/(1024.0f*1024.0f));
  printf("\n avail : %f MB \n",float(avail)/(1024.0f*1024.0f));
  printf("\n used : %f MB \n",float(total-avail)/(1024.0f*1024.0f));

// initialise things if we are looking for recurrences
  if(*RecFLAG==1){
    // create index array for taking norms
    for(int iz=0; iz<2*normWidth; iz++){
      for(int iy=0; iy<2*normWidth; iy++){
	for(int ix=0; ix<normWidth; ix++){
	  int ik = ((*KTZ-normWidth+iz)*(*IKTY)+*KTZ-normWidth+iy)*(*IKTX)+ix;
	  iNorm[(iz*2*normWidth+iy)*normWidth+ix]=ik;
	}
      }
    }
    char Ufile[sizeof "JOB#/UPO_Zk.out"];
    sprintf(Ufile, "JOB%1d/UPO_ZK.out",rank);
    UPOZK = fopen(Ufile,"wb");

    char infofile[sizeof "JOB#/UPOinfo_ts.dat"];
    sprintf(infofile, "JOB%1d/UPOinfo_ts.dat",rank);
    UPOfile = fopen(infofile,"w");

    //        Write header for UPO Guesses output file
    fprintf(UPOfile, "Guess No : Start time : guess Period : guess Shift x: guess Shift y: guess Residual \n");

    // set up shifts
    for(int iShiftY=0; iShiftY < NSHIFTY; iShiftY++){
      shiftY[iShiftY] = twopi*iShiftY/(NSHIFTY);
    }
    for(int iShiftX=0; iShiftX < NSHIFTX; iShiftX++){
      shiftX[iShiftX] = twopi*iShiftX/(NSHIFTX*(*alpha));
    }
    for(int iShiftZ=0; iShiftZ < NSHIFTZ; iShiftZ++){
      shiftZ[iShiftZ] = twopi*iShiftZ/(NSHIFTZ);
    }
    // get norm of I.C.
    normZ = computeNorm(XIIstore,ETAstore,ZETstore,RHOstore,LC,0);
    normZStore[0]=normZ;
  }
