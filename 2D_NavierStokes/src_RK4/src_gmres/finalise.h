/////////////////////////////////////////////////////////////////////////////////////
//  wrap up and write means, then free allocated memory
/////////////////////////////////////////////////////////////////////////////////////
  // Free GPU global memory
  (cudaFree(d_Z));
  (cudaFree(d_Z0));
  (cudaFree(d_Z1));
  (cudaFree(d_Z2));
  (cudaFree(d_Z3));
  (cudaFree(d_UK));
  (cudaFree(d_VK));
  (cudaFree(d_NZK));
  
  (cudaFree(d_ZR));	
  (cudaFree(d_UR));
  (cudaFree(d_VR));
  (cudaFree(d_NZR));
  
  (cudaFree(d_ikF));
  (cudaFree(d_ikN));
  (cudaFree(d_kx));
  (cudaFree(d_ky));
  (cudaFree(d_LL));
  
  //Destroy fft plans
  (cufftDestroy(PlanD2Z));
  (cufftDestroy(PlanZ2D));
  (cufftDestroy(PlanBatchD2Z));


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
