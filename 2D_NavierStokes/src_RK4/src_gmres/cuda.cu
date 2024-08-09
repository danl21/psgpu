#include <stdio.h>
#include <stdlib.h>
#include <cufft.h>
#include <cuda_runtime.h>
#include <thrust/extrema.h>
#include <thrust/reduce.h>
#include <thrust/device_vector.h>
#include "../kernels.h"
#include "../functions.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////

extern "C" void timestep_cuda_(double2 *Z,
			       double *ZR,
			       double *kx,
			       double *ky,
			       double *Time,
			       double *TSTART,
			       double *TNEWT,
			       double *DELT,
			       double *alpha,
			       double *V2,
			       int *ikF,
			       int *ikN,
			       int * IKTX,
			       int * IKTY,
			       int * KTY,
			       int * NX,
			       int * NY,
			       int *NSTOP,
			       int *LC)
			      
{
  // Define CPU variables
  printf("\n In timestep_cuda \n");
  //Place parameters into global holders
  #include "init.h"
  //  printf("\n nstop= %d,  dt= : %e \n",*NSTOP,*DELT);
  
  // Do an initial set of velocity coeffs. Subsequently this occurs at the end of stepping kernels
  setVelocity<<<nblocks,nthreads>>>(d_Z,d_Z0, d_UK, d_VK,d_kx,d_ky,d_LL);

  // **************************
  //STEPPING STARTS HERE
  // **************************
  *Time = 0.0;
  for(int NT=0; NT<*NSTOP; NT++){
    
    KR_FFT_ALL(d_Z0,d_UK,d_VK,d_ZR,d_UR,d_VR,PlanZ2D,d_ikF);
     
    multReal<<<nblocks,nthreads>>>(d_ZR,d_UR,d_VR,d_NZR);    // Do real space convolution, both terms inside NZR 
    RK_FFT(d_NZK,d_NZR,PlanBatchD2Z,d_ikN);    // RK does a batch of 2 ffts, previously NZK and UK
    Step1<<<nblocks,nthreads>>>(d_Z,d_Z0,d_Z1,d_UK,d_VK,d_NZK,d_kx,d_ky,d_LL);
    KR_FFT_ALL(d_Z0,d_UK,d_VK,d_ZR,d_UR,d_VR,PlanZ2D,d_ikF);    // Second half of time step:

    multReal<<<nblocks,nthreads>>>(d_ZR,d_UR,d_VR,d_NZR);
    RK_FFT(d_NZK,d_NZR,PlanBatchD2Z,d_ikN);
    Step1<<<nblocks,nthreads>>>(d_Z,d_Z0,d_Z2,d_UK,d_VK,d_NZK,d_kx,d_ky,d_LL);
    KR_FFT_ALL(d_Z0,d_UK,d_VK,d_ZR,d_UR,d_VR,PlanZ2D,d_ikF);    // Second half of time step:

    multReal<<<nblocks,nthreads>>>(d_ZR,d_UR,d_VR,d_NZR);
    RK_FFT(d_NZK,d_NZR,PlanBatchD2Z,d_ikN);
    Step2<<<nblocks,nthreads>>>(d_Z,d_Z0,d_Z3,d_UK,d_VK,d_NZK,d_kx,d_ky,d_LL);
    KR_FFT_ALL(d_Z0,d_UK,d_VK,d_ZR,d_UR,d_VR,PlanZ2D,d_ikF);    // Second half of time step:

    multReal<<<nblocks,nthreads>>>(d_ZR,d_UR,d_VR,d_NZR);
    RK_FFT(d_NZK,d_NZR,PlanBatchD2Z,d_ikN);
    Step<<<nblocks,nthreads>>>(d_Z,d_Z0,d_Z1,d_Z2,d_Z3,d_UK,d_VK,d_NZK,d_kx,d_ky,d_LL);


    *Time += delt;
    tme = *Time;// Increase time
  }

  // Copy final state off GPU

  (cudaMemcpy(Z, d_Z, sizeof(cufftDoubleComplex)*nkt, cudaMemcpyDeviceToHost));

  #include "finalise.h"
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
