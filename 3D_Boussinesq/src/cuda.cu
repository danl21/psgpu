////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
//    THIS IS THE MAIN CUDA CALLING FUNCTION FOR THE TIMESTEPPING, CALLED BY THE FORTRAN
//    MAIN PROGRAM. THIS IS MAINLY DUE TO THE GMRES CODE BEING WRITTEN IN FORTRAN AND NOT
//    PORTED TO C (YET?). BECAUSE OF THIS, EVERYTHING NEEDS TO BE PASSED TO TIMESTEP_CUDA()
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <time.h>
#include <cufft.h>
#include <helper_cuda.h>
#include <cuda_runtime.h>
#include <thrust/extrema.h>
#include <thrust/device_vector.h>
#include "kernels.h"
#include "functions.h"

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

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
  
  printf("\n In timestep_cuda \n");

#include "init.h" // Set everything up: allocate memory, set various variables for diagnostics etc.

  setVelocity<<<nblocks,nthreads>>>(d_xii,d_eta,d_zet,d_rho,d_x0,d_e0,d_z0,d_r0,d_UK,d_VK,d_WK,d_kx,d_ky,d_kz,d_LL); // get initial velocity given vorticity initial condition
 
  // **************************
  //STEPPING STARTS HERE
  // **************************
  fflush(stdout);
  int count = 0;
  countStats = 0;
  int timestep=0;
  while(tme<=*Tstep){
    timestep++;

    KR_FFT_ALL(d_x0,d_e0,d_z0,d_r0,d_UK,d_VK,d_WK,d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,d_rx,d_ry,d_rz,PlanZ2D,d_ikF,n2); // convert to physical space
    
#include "out_low.h" // low frequency outputs/diagnostics

    multReal<<<nblocks,nthreads>>>(d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,d_rx,d_ry,d_rz); // real space convolution, result stored in U,V,W	 
    RK_FFT(d_UK,d_VK,d_WK,d_rk,d_UR,d_VR,d_WR,d_rx,PlanD2Z,d_ikN,d_LL,n2); // convert to Fourier
    Step1<<<nblocks,nthreads>>>(d_xii,d_eta,d_zet,d_rho,d_x0,d_e0,d_z0,d_r0,d_x1,d_e1,d_z1,d_r1,d_UK,d_VK,d_WK,d_rk,d_kx,d_ky,d_kz,d_LL); // First RK4 substep
    
    KR_FFT_ALL(d_x0,d_e0,d_z0,d_r0,d_UK,d_VK,d_WK,d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,d_rx,d_ry,d_rz,PlanZ2D,d_ikF,n2);// convert to physical space
    multReal<<<nblocks,nthreads>>>(d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,d_rx,d_ry,d_rz); // real space convolution, result stored in U,V,W
    RK_FFT(d_UK,d_VK,d_WK,d_rk,d_UR,d_VR,d_WR,d_rx,PlanD2Z,d_ikN,d_LL,n2); // convert to Fourier
    Step1<<<nblocks,nthreads>>>(d_xii,d_eta,d_zet,d_rho,d_x0,d_e0,d_z0,d_r0,d_x2,d_e2,d_z2,d_r2,d_UK,d_VK,d_WK,d_rk,d_kx,d_ky,d_kz,d_LL); // Second RK4 substep

    KR_FFT_ALL(d_x0,d_e0,d_z0,d_r0,d_UK,d_VK,d_WK,d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,d_rx,d_ry,d_rz,PlanZ2D,d_ikF,n2);// convert to physical space
    multReal<<<nblocks,nthreads>>>(d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,d_rx,d_ry,d_rz); // real space convolution, result stored in U,V,W	 
    RK_FFT(d_UK,d_VK,d_WK,d_rk,d_UR,d_VR,d_WR,d_rx,PlanD2Z,d_ikN,d_LL,n2);// convert to Fourier
    Step2<<<nblocks,nthreads>>>(d_xii,d_eta,d_zet,d_rho,d_x0,d_e0,d_z0,d_r0,d_x3,d_e3,d_z3,d_r3,d_UK,d_VK,d_WK,d_rk,d_kx,d_ky,d_kz,d_LL);// Third RK4 substep

    KR_FFT_ALL(d_x0,d_e0,d_z0,d_r0,d_UK,d_VK,d_WK,d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,d_rx,d_ry,d_rz,PlanZ2D,d_ikF,n2);// convert to physical space
    multReal<<<nblocks,nthreads>>>(d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,d_rx,d_ry,d_rz); // real space convolution, result stored in U,V,W	 
    RK_FFT(d_UK,d_VK,d_WK,d_rk,d_UR,d_VR,d_WR,d_rx,PlanD2Z,d_ikN,d_LL,n2);// convert to Fourier
    Step<<<nblocks,nthreads>>>(d_xii,d_eta,d_zet,d_rho,d_x0,d_e0,d_z0,d_r0,d_x1,d_e1,d_z1,d_r1,d_x2,d_e2,d_z2,d_r2,d_x3,d_e3,d_z3,d_r3,d_UK,d_VK,d_WK,d_rk,d_kx,d_ky,d_kz,d_LL);// Final RK4 step

    *Time = *Time + delt;    // Increase time
    tme += delt;
    
#include "out_high.h" // High frequency outputs/diagnostics
    if(*Dtarget != 0){
      (cudaMemcpyToSymbol(d_AMPFOR,&ampfor,sizeof(double))); // if we are throttling copy the forcing amplitude to the device
    }
  }
  fflush(stdout);

    // Copy state off GPU
  (cudaMemcpy(xii, d_xii, sizeof(cufftDoubleComplex)*nkt, cudaMemcpyDeviceToHost));
  (cudaMemcpy(eta, d_eta, sizeof(cufftDoubleComplex)*nkt, cudaMemcpyDeviceToHost));
  (cudaMemcpy(zet, d_zet, sizeof(cufftDoubleComplex)*nkt, cudaMemcpyDeviceToHost));
  (cudaMemcpy(rho, d_rho, sizeof(cufftDoubleComplex)*nkt, cudaMemcpyDeviceToHost));

#include "finalise.h" // wrap up: do final diagnostics (pdfs, means, rmses etc) and free all allocated memory
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
