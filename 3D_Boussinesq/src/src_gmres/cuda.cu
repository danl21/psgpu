#include <stdio.h>
#include <time.h>
#include <cufft.h>
#include <helper_cuda.h>
#include <cuda_runtime.h>
#include "../kernels.h"
#include "../functions.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\
////////////////////////

extern "C" void timestep_cuda_(double2 *xii,
			       double2 *eta,
			       double2 *zet,
			       double2 *rho,
			       double *kx,
			       double *ky,
			       double *kz,
			       double *Time, 
			       double *TSTART,
			       double *Tstep,
			       double *AMPFOR,
			       double *DELT,
			       double *ResThresh,
			       double *KFY,
			       double *V2,
			       double *Ri,
			       double *Sc,
			       double *Theta,
			       double *alpha,
			       double *Dtarget,
			       int *ikF,
			       int *ikN,
			       int *LC,
			       int * IKTX,
			       int * IKTY,
			       int * IKTZ,
			       int * KTZ,
			       int *NKT,
			       int * NX,
			       int * NY,
			       int *NZ,
			       int *NOUT,
			       int *NOUTV,
			       int *statsFLAG,
			       int *RecFLAG,
			       int *adaptFLAG,
			       int *RANK
			       )
{
  
  // Define CPU variables
  printf("\n In timestep_cuda \n");
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
  delt = *DELT;
  nOut = *NOUT;
  tme = *TSTART;
  tstart = *TSTART;
  ampfor = -*AMPFOR;
  v2 = *V2;
  ResidualThreshold = *ResThresh;
  rank = *RANK;
  printf("Tstep %e \n",*Tstep);
  printf("ampfor %e \n",ampfor);
  printf("delt %e \n",delt);
  printf("v2 %e \n",v2);
  printf("Ri %e\n",*Ri);
  printf("Sc %e\n",*Sc);
  printf("theta %e \n",*Theta);
  printf("kfy %e \n",kfy);

  double IN = 1.0/nr;
  sinTh = sin(*Theta);
  cosTh = cos(*Theta);
  size_t avail,total;
  
  // Define global device variables
  int *d_ikF, *d_ikN, *d_LL;
  cufftDoubleComplex *d_xii,*d_eta,*d_zet,*d_rho,*d_UK, *d_VK, *d_WK,*d_rk;
  cufftDoubleComplex *d_x0,*d_x1,*d_x2,*d_x3,*d_e0,*d_e1,*d_e2,*d_e3,*d_z0,*d_z1,*d_z2,*d_z3,*d_r0,*d_r1,*d_r2,*d_r3;
  cufftDoubleReal *d_xir,*d_etr,*d_ztr,*d_rx,*d_ry,*d_rz,*d_UR,*d_VR,*d_WR,*d_ro;
  cufftHandle PlanZ2D,PlanD2Z;

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
  // Do a check of global memory use
  avail =0;
  total = 0;
  cudaMemGetInfo(&avail,&total);
  
  // Do an initial set of velocity coeffs. Subsequently this occurs at the end of stepping kernels
  setVelocity<<<nblocks,nthreads>>>(d_xii,d_eta,d_zet,d_rho,d_x0,d_e0,d_z0,d_r0,d_UK,d_VK,d_WK,d_kx,d_ky,d_kz,d_LL);
 
  // **************************
  //STEPPING STARTS HERE
  // **************************
  fflush(stdout);
  int timestep=0;
  while(*Tstep-tme > delt){
    timestep++;

    KR_FFT_ALL(d_x0,d_e0,d_z0,d_r0,d_UK,d_VK,d_WK,d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,d_rx,d_ry,d_rz,PlanZ2D,d_ikF,n2);

    multReal<<<nblocks,nthreads>>>(d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,d_rx,d_ry,d_rz); // real space convolution, result in U,V,W	 
    RK_FFT(d_UK,d_VK,d_WK,d_rk,d_UR,d_VR,d_WR,d_rx,PlanD2Z,d_ikN,d_LL,n2);
    Step1<<<nblocks,nthreads>>>(d_xii,d_eta,d_zet,d_rho,d_x0,d_e0,d_z0,d_r0,d_x1,d_e1,d_z1,d_r1,d_UK,d_VK,d_WK,d_rk,d_kx,d_ky,d_kz,d_LL);
    
    KR_FFT_ALL(d_x0,d_e0,d_z0,d_r0,d_UK,d_VK,d_WK,d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,d_rx,d_ry,d_rz,PlanZ2D,d_ikF,n2);
    multReal<<<nblocks,nthreads>>>(d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,d_rx,d_ry,d_rz); // real space convolution, result in U,V,W
    RK_FFT(d_UK,d_VK,d_WK,d_rk,d_UR,d_VR,d_WR,d_rx,PlanD2Z,d_ikN,d_LL,n2);
    Step1<<<nblocks,nthreads>>>(d_xii,d_eta,d_zet,d_rho,d_x0,d_e0,d_z0,d_r0,d_x2,d_e2,d_z2,d_r2,d_UK,d_VK,d_WK,d_rk,d_kx,d_ky,d_kz,d_LL);	 

    KR_FFT_ALL(d_x0,d_e0,d_z0,d_r0,d_UK,d_VK,d_WK,d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,d_rx,d_ry,d_rz,PlanZ2D,d_ikF,n2);
    multReal<<<nblocks,nthreads>>>(d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,d_rx,d_ry,d_rz); // real space convolution, result in U,V,W
    RK_FFT(d_UK,d_VK,d_WK,d_rk,d_UR,d_VR,d_WR,d_rx,PlanD2Z,d_ikN,d_LL,n2);
    Step2<<<nblocks,nthreads>>>(d_xii,d_eta,d_zet,d_rho,d_x0,d_e0,d_z0,d_r0,d_x3,d_e3,d_z3,d_r3,d_UK,d_VK,d_WK,d_rk,d_kx,d_ky,d_kz,d_LL);	 

    KR_FFT_ALL(d_x0,d_e0,d_z0,d_r0,d_UK,d_VK,d_WK,d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,d_rx,d_ry,d_rz,PlanZ2D,d_ikF,n2);
    multReal<<<nblocks,nthreads>>>(d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,d_rx,d_ry,d_rz); // real space convolution, result in U,V,W
    RK_FFT(d_UK,d_VK,d_WK,d_rk,d_UR,d_VR,d_WR,d_rx,PlanD2Z,d_ikN,d_LL,n2);
    Step<<<nblocks,nthreads>>>(d_xii,d_eta,d_zet,d_rho,d_x0,d_e0,d_z0,d_r0,d_x1,d_e1,d_z1,d_r1,d_x2,d_e2,d_z2,d_r2,d_x3,d_e3,d_z3,d_r3,d_UK,d_VK,d_WK,d_rk,d_kx,d_ky,d_kz,d_LL);	 

    *Time = *Time + delt;    // Increase time
    tme += delt;
    
  }
  if( tme != *Tstep){
    delt = *Tstep-tme;
    printf("adjusting final timestep dt = %e \n",delt);
    (cudaMemcpyToSymbol(d_DELT,&delt,sizeof(double)));
    KR_FFT_ALL(d_x0,d_e0,d_z0,d_r0,d_UK,d_VK,d_WK,d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,d_rx,d_ry,d_rz,PlanZ2D,d_ikF,n2);

    multReal<<<nblocks,nthreads>>>(d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,d_rx,d_ry,d_rz); // real space convolution, result in U,V,W	 
    RK_FFT(d_UK,d_VK,d_WK,d_rk,d_UR,d_VR,d_WR,d_rx,PlanD2Z,d_ikN,d_LL,n2);
    Step1<<<nblocks,nthreads>>>(d_xii,d_eta,d_zet,d_rho,d_x0,d_e0,d_z0,d_r0,d_x1,d_e1,d_z1,d_r1,d_UK,d_VK,d_WK,d_rk,d_kx,d_ky,d_kz,d_LL);
    
    KR_FFT_ALL(d_x0,d_e0,d_z0,d_r0,d_UK,d_VK,d_WK,d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,d_rx,d_ry,d_rz,PlanZ2D,d_ikF,n2);
    multReal<<<nblocks,nthreads>>>(d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,d_rx,d_ry,d_rz); // real space convolution, result in U,V,W
    RK_FFT(d_UK,d_VK,d_WK,d_rk,d_UR,d_VR,d_WR,d_rx,PlanD2Z,d_ikN,d_LL,n2);
    Step1<<<nblocks,nthreads>>>(d_xii,d_eta,d_zet,d_rho,d_x0,d_e0,d_z0,d_r0,d_x2,d_e2,d_z2,d_r2,d_UK,d_VK,d_WK,d_rk,d_kx,d_ky,d_kz,d_LL);	 

    KR_FFT_ALL(d_x0,d_e0,d_z0,d_r0,d_UK,d_VK,d_WK,d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,d_rx,d_ry,d_rz,PlanZ2D,d_ikF,n2);
    multReal<<<nblocks,nthreads>>>(d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,d_rx,d_ry,d_rz); // real space convolution, result in U,V,W
    RK_FFT(d_UK,d_VK,d_WK,d_rk,d_UR,d_VR,d_WR,d_rx,PlanD2Z,d_ikN,d_LL,n2);
    Step2<<<nblocks,nthreads>>>(d_xii,d_eta,d_zet,d_rho,d_x0,d_e0,d_z0,d_r0,d_x3,d_e3,d_z3,d_r3,d_UK,d_VK,d_WK,d_rk,d_kx,d_ky,d_kz,d_LL);	 

    KR_FFT_ALL(d_x0,d_e0,d_z0,d_r0,d_UK,d_VK,d_WK,d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,d_rx,d_ry,d_rz,PlanZ2D,d_ikF,n2);
    multReal<<<nblocks,nthreads>>>(d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,d_rx,d_ry,d_rz); // real space convolution, result in U,V,W
    RK_FFT(d_UK,d_VK,d_WK,d_rk,d_UR,d_VR,d_WR,d_rx,PlanD2Z,d_ikN,d_LL,n2);
    Step<<<nblocks,nthreads>>>(d_xii,d_eta,d_zet,d_rho,d_x0,d_e0,d_z0,d_r0,d_x1,d_e1,d_z1,d_r1,d_x2,d_e2,d_z2,d_r2,d_x3,d_e3,d_z3,d_r3,d_UK,d_VK,d_WK,d_rk,d_kx,d_ky,d_kz,d_LL);	 

  }
  fflush(stdout);

    // Copy state off GPU
  (cudaMemcpy(xii, d_xii, sizeof(cufftDoubleComplex)*nkt, cudaMemcpyDeviceToHost));
  (cudaMemcpy(eta, d_eta, sizeof(cufftDoubleComplex)*nkt, cudaMemcpyDeviceToHost));
  (cudaMemcpy(zet, d_zet, sizeof(cufftDoubleComplex)*nkt, cudaMemcpyDeviceToHost));
  (cudaMemcpy(rho, d_rho, sizeof(cufftDoubleComplex)*nkt, cudaMemcpyDeviceToHost));

  printf("pdfs and averages done \n"); 
  // Free GPU global memory
  (cudaFree(d_xii));
  (cudaFree(d_eta));
  (cudaFree(d_zet));
  (cudaFree(d_rho));
  (cudaFree(d_UK));
  (cudaFree(d_VK));
  (cudaFree(d_WK));
  (cudaFree(d_rk));
  
  (cudaFree(d_x0));	
  (cudaFree(d_e0));	
  (cudaFree(d_z0));	
  (cudaFree(d_r0));	

  (cudaFree(d_x1));	
  (cudaFree(d_e1));	
  (cudaFree(d_z1));	
  (cudaFree(d_r1));	

  (cudaFree(d_x2));	
  (cudaFree(d_e2));	
  (cudaFree(d_z2));	
  (cudaFree(d_r2));	

  (cudaFree(d_x3));	
  (cudaFree(d_e3));	
  (cudaFree(d_z3));	
  (cudaFree(d_r3));	

  (cudaFree(d_xir));	
  (cudaFree(d_etr));	
  (cudaFree(d_ztr));	
  (cudaFree(d_rx));	
  (cudaFree(d_ry));	
  (cudaFree(d_rz));	
  (cudaFree(d_ro));	
  (cudaFree(d_UR));
  (cudaFree(d_VR));
  (cudaFree(d_WR));
  
  (cudaFree(d_LL));
  (cudaFree(d_ikF));
  (cudaFree(d_ikN));
  (cudaFree(d_kx));
  (cudaFree(d_ky));
  (cudaFree(d_kz));
  printf("GPU free done \n"); 
  //Destroy fft plans
  (cufftDestroy(PlanZ2D));
  (cufftDestroy(PlanD2Z));
   
  printf("time stepping done \n"); 
  fflush(stdout);

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\
////////////////////////

extern "C" void timederiv_cuda_(double2 *xii,
			       double2 *eta,
			       double2 *zet,
			       double2 *rho,
			       double *kx,
			       double *ky,
			       double *kz,
			       double *Time, 
			       double *TSTART,
			       double *Tstep,
			       double *AMPFOR,
			       double *DELT,
			       double *ResThresh,
			       double *KFY,
			       double *V2,
			       double *Ri,
			       double *Sc,
			       double *Theta,
			       double *alpha,
			       double *Dtarget,
			       int *ikF,
			       int *ikN,
			       int *LC,
			       int * IKTX,
			       int * IKTY,
			       int * IKTZ,
			       int * KTZ,
			       int *NKT,
			       int * NX,
			       int * NY,
			       int *NZ,
			       int *NOUT,
			       int *NOUTV,
			       int *statsFLAG,
			       int *RecFLAG,
			       int *adaptFLAG,
			       int *RANK
			       )
{
  
  // Define CPU variables
  printf("\n In timestep_cuda \n");
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
  double nStop = max(*Tstep/(*DELT),1.0);
  if(nStop**DELT != *Tstep){
    delt = *Tstep/nStop;
      }else{
    delt = *DELT;  
  }
  nOut = *NOUT;
  tme = *TSTART;
  tstart = *TSTART;
  ampfor = -*AMPFOR;
  v2 = *V2;
  ResidualThreshold = *ResThresh;
  rank = *RANK;
  printf("Tstep %e \n",*Tstep);
  printf("ampfor %e \n",ampfor);
  printf("delt %e \n",delt);
  printf("v2 %e \n",v2);
  printf("Ri %e\n",*Ri);
  printf("Sc %e\n",*Sc);
  printf("theta %e \n",*Theta);
  printf("kfy %e \n",kfy);

  double IN = 1.0/nr;
  sinTh = sin(*Theta);
  cosTh = cos(*Theta);
  printf("sin(theta) %e \n",sinTh);
  printf("cos(theta) %e \n",cosTh);
  
  // Define global device variables
  int *d_ikF, *d_ikN, *d_LL;
  cufftDoubleComplex *d_xii,*d_eta,*d_zet,*d_rho,*d_UK, *d_VK, *d_WK,*d_rk;
  cufftDoubleComplex *d_x0,*d_e0,*d_z0,*d_r0;
  cufftDoubleReal *d_xir,*d_etr,*d_ztr,*d_rx,*d_ry,*d_rz,*d_UR,*d_VR,*d_WR,*d_ro;
  cufftHandle PlanZ2D,PlanD2Z;

  // Allocate global memory on GPU. (Constant memory does not need allocating) 	
  (cudaMalloc((void**)&d_xii,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_x0,sizeof(cufftDoubleComplex)*(nkt)));

  (cudaMalloc((void**)&d_eta,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_e0,sizeof(cufftDoubleComplex)*(nkt)));

  (cudaMalloc((void**)&d_zet,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_z0,sizeof(cufftDoubleComplex)*(nkt)));

  (cudaMalloc((void**)&d_rho,sizeof(cufftDoubleComplex)*(nkt)));
  (cudaMalloc((void**)&d_r0,sizeof(cufftDoubleComplex)*(nkt)));

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
  // Do an initial set of velocity coeffs. Subsequently this occurs at the end of stepping kernels
  setVelocity<<<nblocks,nthreads>>>(d_xii,d_eta,d_zet,d_rho,d_x0,d_e0,d_z0,d_r0,d_UK,d_VK,d_WK,d_kx,d_ky,d_kz,d_LL);
  KR_FFT_ALL(d_x0,d_e0,d_z0,d_r0,d_UK,d_VK,d_WK,d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,d_rx,d_ry,d_rz,PlanZ2D,d_ikF,n2);
  
  multReal<<<nblocks,nthreads>>>(d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,d_rx,d_ry,d_rz); // real space convolution, result in U,V,W	 
  RK_FFT(d_UK,d_VK,d_WK,d_rk,d_UR,d_VR,d_WR,d_rx,PlanD2Z,d_ikN,d_LL,n2);
  RHS<<<nblocks,nthreads>>>(d_xii,d_eta,d_zet,d_rho,d_x0,d_e0,d_z0,d_r0,d_UK,d_VK,d_WK,d_rk,d_kx,d_ky,d_kz,d_LL);
  cudaThreadSynchronize();
    // Copy state off GPU
  (cudaMemcpy(xii, d_x0, sizeof(cufftDoubleComplex)*nkt, cudaMemcpyDeviceToHost));
  (cudaMemcpy(eta, d_e0, sizeof(cufftDoubleComplex)*nkt, cudaMemcpyDeviceToHost));
  (cudaMemcpy(zet, d_z0, sizeof(cufftDoubleComplex)*nkt, cudaMemcpyDeviceToHost));
  (cudaMemcpy(rho, d_r0, sizeof(cufftDoubleComplex)*nkt, cudaMemcpyDeviceToHost));

  printf("pdfs and averages done \n"); 
  // Free GPU global memory
  (cudaFree(d_xii));
  (cudaFree(d_eta));
  (cudaFree(d_zet));
  (cudaFree(d_rho));
  (cudaFree(d_UK));
  (cudaFree(d_VK));
  (cudaFree(d_WK));
  (cudaFree(d_rk));
  
  (cudaFree(d_x0));	
  (cudaFree(d_e0));	
  (cudaFree(d_z0));	
  (cudaFree(d_r0));	

  (cudaFree(d_xir));	
  (cudaFree(d_etr));	
  (cudaFree(d_ztr));	
  (cudaFree(d_rx));	
  (cudaFree(d_ry));	
  (cudaFree(d_rz));	
  (cudaFree(d_ro));	
  (cudaFree(d_UR));
  (cudaFree(d_VR));
  (cudaFree(d_WR));
  
  (cudaFree(d_LL));
  (cudaFree(d_ikF));
  (cudaFree(d_ikN));
  (cudaFree(d_kx));
  (cudaFree(d_ky));
  (cudaFree(d_kz));
  //Destroy fft plans
  (cufftDestroy(PlanZ2D));
  (cufftDestroy(PlanD2Z));
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
