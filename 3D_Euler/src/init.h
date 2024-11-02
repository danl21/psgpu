  int nVort = *NOUTV;
  nOut = *NOUT;
  ny = *NY;
  nx = *NX;
  nz = *NZ;
  nkt = *NKT; 
  nr  = (*NY)*(*NX)*(*NZ);
  nx2 = (*NX)/2 +1;
  n2 = nx2*(*NY)*(*NZ);
tme = *TSTART;
tstart = *TSTART;

  freq2 = 2*nOut*delt;
  tavgCount =0;
  double2 *XIItavg,*ETAtavg,*ZETtavg;
  double  *Umean,*Vmean,*Wmean,*uhat,*vhat,*what;

  double IN = 1.0/nr;
  size_t avail,total;
  
  // Define global device variables
  int *d_ikF, *d_ikN, *d_LL;
  cufftDoubleComplex *d_xii,*d_eta,*d_zet, *d_UK, *d_VK, *d_WK, *d_RHX,*d_RHY,*d_RHZ;
  cufftDoubleReal *d_xir,*d_etr,*d_ztr, *d_UR,*d_VR,*d_WR;
  double *d_Urms,*d_Vrms,*d_Wrms;
  cufftHandle PlanZ2D,PlanD2Z;

  // Allocate global memory on GPU. (Constant memory does not need allocating) 	
  (cudaMalloc((void**)&d_xii,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_eta,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_zet,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_UK,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_VK,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_WK,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_RHX,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_RHY,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_RHZ,sizeof(cufftDoubleComplex)*(nkt)));
  
  (cudaMalloc((void**)&d_xir,sizeof(cufftDoubleReal)*(nr)));
  (cudaMalloc((void**)&d_etr,sizeof(cufftDoubleReal)*(nr)));
  (cudaMalloc((void**)&d_ztr,sizeof(cufftDoubleReal)*(nr)));
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
  (cudaMemcpyToSymbol(d_DELT,DELT,sizeof(double)));
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
  printf("NX= %d \n",*NX);
  printf("NZ= %d \n",*NZ);
  fflush(stdout);

  char xiifile[sizeof "xii000.raw"];
  char etafile[sizeof "eta000.raw"];
  char zetfile[sizeof "zet000.raw"];

  // Set up diagnostic files and arrays
  int nStats = *NSTOP/nOut+1;
  printf("nStats= %d \n",nStats);
  fflush(stdout);
  if(*statsFLAG !=0){
    char statsfile[sizeof "stats.dat"];
    sprintf(statsfile, "stats.dat");
    stats = fopen(statsfile,"w");

    char engyfile[sizeof "engy.dat"];
    sprintf(engyfile, "engy.dat");
    engy = fopen(engyfile,"w");

    char avgsfile[sizeof "avgs.dat"];
    sprintf(avgsfile, "avgs.dat");
    avgs = fopen(avgsfile,"w");

    char meansfile[sizeof "means.dat"];
    sprintf(meansfile, "means.dat");
    means = fopen(meansfile,"w");

    char pdfsfile[sizeof "pdfs.dat"];
    sprintf(pdfsfile, "pdfs.dat");
    pdfs = fopen(pdfsfile,"w");

    XIItavg=(double2*)malloc(sizeof(double2)*nkt);
    ETAtavg=(double2*)malloc(sizeof(double2)*nkt);
    ZETtavg=(double2*)malloc(sizeof(double2)*nkt);

    (cudaMalloc((void**)&d_Urms,sizeof(cufftDoubleReal)*nr));
    (cudaMalloc((void**)&d_Vrms,sizeof(cufftDoubleReal)*nr));
    (cudaMalloc((void**)&d_Wrms,sizeof(cufftDoubleReal)*nr));

    Energy=(double*)malloc(sizeof(double)*nStats);

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
  }
  // Do a check of global memory use
  avail =0;
  total = 0;
  cudaMemGetInfo(&avail,&total);
  
  printf("\n total : %f MB \n",float(total)/(1024.0f*1024.0f));
  printf("\n avail : %f MB \n",float(avail)/(1024.0f*1024.0f));
  printf("\n used : %f MB \n",float(total-avail)/(1024.0f*1024.0f));


