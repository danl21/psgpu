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
#include <cuda_runtime.h>
#include "kernels.h"
#include "functions.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

extern "C" void timestep_cuda_(double2 *xii,
			       double2 *eta,
			       double2 *zet,
			       double *xir,
			       double *etr,
			       double *ztr,
			       double *kx,
			       double *ky,
			       double *kz,
			       double *Time, 
			       double *TSTART,
			       double *DELT,
			       double *alpha,
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
			       int *NSTOP,
			       int *NOUT,
			       int *NOUTV,
			       int *statsFLAG){

#include "init.h"  // Set everything up: allocate memory, set various variables for diagnostics etc.

  // Do an initial set of velocity coeffs. Subsequently this occurs at the end of stepping kernels
  setVelocity<<<nblocks,nthreads>>>(d_xii,d_eta,d_zet,d_UK,d_VK,d_WK,d_kx,d_ky,d_kz,d_LL);
 
  // **************************
  //STEPPING STARTS HERE
  // **************************
  fflush(stdout);
  int count = 0;
  for(NT=0; NT<*NSTOP; NT++){
    KR_FFT_ALL(d_xii,d_eta,d_zet,d_UK,d_VK,d_WK,d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,PlanZ2D,d_ikF,n2);

    if((NT % (nOut)) == 0 && *statsFLAG !=0) avgVelocity<<<nblocks,nthreads>>>(d_UR,d_VR,d_WR,d_Urms,d_Vrms,d_Wrms,tavgCount);

    if((NT%nVort)==0 && *NSTOP > 1){ // low frequency outputs of the physical space 3D volumes
      count++;
      sprintf(xiifile, "xii%03d.raw",count);
      sprintf(etafile, "eta%03d.raw",count);
      sprintf(zetfile, "zet%03d.raw",count);
      xiiStream = fopen(xiifile,"wb");
      etaStream = fopen(etafile,"wb");
      zetStream = fopen(zetfile,"wb");
      (cudaMemcpy(ztr, d_ztr, sizeof(double)*nr, cudaMemcpyDeviceToHost));
      (cudaMemcpy(xir, d_xir, sizeof(double)*nr, cudaMemcpyDeviceToHost));
      (cudaMemcpy(etr, d_etr, sizeof(double)*nr, cudaMemcpyDeviceToHost));

      fwrite(xir,sizeof(double),nr,xiiStream);
      fwrite(etr,sizeof(double),nr,etaStream);
      fwrite(ztr,sizeof(double),nr,zetStream);
      fclose(etaStream);
      fflush(etaStream);
      fclose(zetStream);
      fflush(zetStream);
      fclose(xiiStream);
      fflush(xiiStream);
    }


    multReal<<<nblocks,nthreads>>>(d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR); // real space convolution, result in U,V,W
	 
    RK_FFT(d_UK,d_VK,d_WK,d_UR,d_VR,d_WR,PlanD2Z,d_ikN,n2);
	  
    // Predictor step: prestep now does the end of 'convol', the step and resets velocity coeffs
    preStep<<<nblocks,nthreads>>>(d_xii,d_eta,d_zet,d_UK,d_VK,d_WK,d_RHX,d_RHY,d_RHZ,d_kx,d_ky,d_kz,d_LL);
    
    KR_FFT_ALL(d_xii,d_eta,d_zet,d_UK,d_VK,d_WK,d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,PlanZ2D,d_ikF,n2);

    multReal<<<nblocks,nthreads>>>(d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR); // real space convolution, result in U,V,W
	 
    RK_FFT(d_UK,d_VK,d_WK,d_UR,d_VR,d_WR,PlanD2Z,d_ikN,n2);
    *Time = *TSTART + (NT+1)*(*DELT);    // Increase time
    tme += *DELT;

    //Correction step: corStep is analagous in structure to prestep.
    corStep<<<nblocks,nthreads>>>(d_xii,d_eta,d_zet,d_UK,d_VK,d_WK,d_RHX,d_RHY,d_RHZ,d_kx,d_ky,d_kz,d_LL);
    
    if((NT % (nOut)) == 0 && *statsFLAG!=0 && NT!=0){ //High frequency outputs, statistics

      cudaDeviceSynchronize();

      (cudaMemcpy(xii,d_xii,sizeof(cufftDoubleComplex)*nkt,cudaMemcpyDeviceToHost));
      (cudaMemcpy(eta,d_eta,sizeof(cufftDoubleComplex)*nkt,cudaMemcpyDeviceToHost));
      (cudaMemcpy(zet,d_zet,sizeof(cufftDoubleComplex)*nkt,cudaMemcpyDeviceToHost));

      if(*statsFLAG==1)computeStats(xii,eta,zet,XIItavg,ETAtavg,ZETtavg,kx,ky,kz,LC);      
      if(*statsFLAG==2 && tme > 500.0 )computeStats(xii,eta,zet,XIItavg,ETAtavg,ZETtavg,kx,ky,kz,LC);      
      fflush(stats);
    }

  }

#include "finalise.h"
   
  printf("time stepping done \n"); 
  fflush(stdout);

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
