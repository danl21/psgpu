//////////////////////////////////////////////////////////////////////
// Set up for cuda timestep functions
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//Place parameters into global holders
  iktx = *IKTX;
  ikty = *IKTY;
  kty = *KTY;
  ny = *NY;
  nx = *NX;
  nkt = (*IKTX)*(*IKTY);
  nr  = (*NY)*(*NX);
  nx2 = (*NX)/2 +1;
  n2 = nx2*(*NY);
  v2 = *V2;
  in = 1.0/nr;
  tstart = *TSTART;
  tme = 0.0;
  delt = *DELT;
  np[0] = *NY;
  np[1] = *NX; 
 
  // Define global device variables
  // For LL couldn't read in fortran logical to C bool
  int *d_LL;
  int *d_ikF, *d_ikN;
  cufftDoubleComplex *d_Z, *d_UK, *d_VK, *d_NZK;
  cufftDoubleComplex *d_Z0, *d_Z1, *d_Z2, *d_Z3;
  cufftDoubleReal*d_ZR,*d_NZR, *d_UR,*d_VR;
  double *d_tme;
  cufftHandle PlanZ2D,PlanBatchD2Z,PlanD2Z;

  // Allocate global memory on GPU. (Constant memory does not need allocating) 	
  (cudaMalloc((void**)&d_Z,sizeof(cufftDoubleComplex)*nkt));
  (cudaMalloc((void**)&d_Z0,sizeof(cufftDoubleComplex)*nkt));
  (cudaMalloc((void**)&d_Z1,sizeof(cufftDoubleComplex)*nkt));
  (cudaMalloc((void**)&d_Z2,sizeof(cufftDoubleComplex)*nkt));
  (cudaMalloc((void**)&d_Z3,sizeof(cufftDoubleComplex)*nkt));
  (cudaMalloc((void**)&d_UK,sizeof(cufftDoubleComplex)*nkt));
  (cudaMalloc((void**)&d_VK,sizeof(cufftDoubleComplex)*nkt));
  (cudaMalloc((void**)&d_NZK,sizeof(cufftDoubleComplex)*2*nkt));
  
  (cudaMalloc((void**)&d_ZR,sizeof(cufftDoubleReal)*nr));
  (cudaMalloc((void**)&d_UR,sizeof(cufftDoubleReal)*nr));
  (cudaMalloc((void**)&d_VR,sizeof(cufftDoubleReal)*nr));
  (cudaMalloc((void**)&d_NZR,sizeof(cufftDoubleReal)*2*nr));
  
  (cudaMalloc((void**)&d_kx,sizeof(double)*nkt));
  (cudaMalloc((void**)&d_ky,sizeof(double)*nkt));
  (cudaMalloc((void**)&d_ikF,sizeof(int)*n2));
  (cudaMalloc((void**)&d_ikN,sizeof(int)*n2));
  (cudaMalloc((void**)&d_LL,sizeof(int)*nkt));
  (cudaMalloc((void**)&d_tme,sizeof(double)));

// Copy state data to GPU global memory 
 (cudaMemcpy(d_Z,Z,sizeof(cufftDoubleComplex)*nkt,cudaMemcpyHostToDevice));
 // (cudaMemcpy(d_ZR,ZR,sizeof(cufftDoubleReal)*nr,cudaMemcpyHostToDevice));
  
  // Copy constant parameters to GPU constant memory
  (cudaMemcpyToSymbol(d_IN,&in,sizeof(double)));
  (cudaMemcpyToSymbol(d_DELT,&delt,sizeof(double)));
  (cudaMemcpyToSymbol(d_NU,V2,sizeof(double)));
  (cudaMemcpyToSymbol(d_IKTX,IKTX,sizeof(int)));
  (cudaMemcpyToSymbol(d_IKTY,IKTY,sizeof(int)));
  (cudaMemcpyToSymbol(d_KTY,KTY,sizeof(int)));
  (cudaMemcpyToSymbol(d_NX,NX,sizeof(int)));
  (cudaMemcpyToSymbol(d_NX2,&nx2,sizeof(int)));
  (cudaMemcpyToSymbol(d_NY,NY,sizeof(int)));
  (cudaMemcpyToSymbol(d_OR,&nr,sizeof(int)));
  (cudaMemcpyToSymbol(d_OK,&nkt,sizeof(int)));
  (cudaMemcpyToSymbol(d_O2,&n2,sizeof(int)));

  //Set up various arrays to enable generic kernel calls
  //i..e=/ calculate indexing for padding either side of FFTs, wavenumber arrays, mask, and timestep arrays.
  // This must be done on CPU for scalability (large problems violate max threads per block)
  (cudaMemcpy(d_kx,kx,sizeof(double)*nkt,cudaMemcpyHostToDevice));
  (cudaMemcpy(d_ky,ky,sizeof(double)*nkt,cudaMemcpyHostToDevice));
  (cudaMemcpy(d_LL,LC,sizeof(int)*nkt,cudaMemcpyHostToDevice));
  (cudaMemcpy(d_ikF,ikF,sizeof(int)*n2,cudaMemcpyHostToDevice));
  (cudaMemcpy(d_ikN,ikN,sizeof(int)*n2,cudaMemcpyHostToDevice));
  //Set FFT Plans (Batched Z2D fails on large problems, D2Z is ok. At least up to 4096^2)
  (cufftPlan2d(&PlanZ2D,*NY,*NX,CUFFT_Z2D));
  (cufftPlan2d(&PlanD2Z,*NY,*NX,CUFFT_D2Z));
  (cufftPlanMany(&PlanBatchD2Z,2,np,NULL,1,0,NULL,1,0,CUFFT_D2Z,2));

printf("T= %16.12e, DT=%16.12e, nStop=%d \n",*Time, delt,*NSTOP);

