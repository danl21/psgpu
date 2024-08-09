/////////////////////////////////////////////////////////////////////////////////////
//  wrap up and write means, then free allocated memory
/////////////////////////////////////////////////////////////////////////////////////

  //  subLam<<<1,1>>>(d_Z,v2);
  KR_FFT_ALL(d_Z,d_UK,d_VK,d_ZR,d_UR,d_VR,PlanZ2D,d_ikF);
  (cudaMemcpy(ZR,(double *) d_ZR, sizeof(double)*nr, cudaMemcpyDeviceToHost));
  fwrite(Time,sizeof(double),1,vort);
  fwrite(ZR,sizeof(double),nr,vort);
  fflush(vort);

  if(*statsFLAG==1){
    (cudaMemcpy(d_Z,ZKtavg,sizeof(cufftDoubleComplex)*nkt,cudaMemcpyHostToDevice));  
    setVelocity<<<nblocks,nthreads>>>(d_Z,d_Z0, d_UK, d_VK,d_kx,d_ky,d_LL);
    KR_FFT_ALL(d_Z,d_UK,d_VK,d_ZR,d_UR,d_VR,PlanZ2D,d_ikF);
    
    (cudaMemcpy(ZR,(double *) d_UR, sizeof(double)*nr, cudaMemcpyDeviceToHost));
    fwrite(Time,sizeof(double),1,avgs);
    fwrite(ZR,sizeof(double),nr,avgs);
    
    (cudaMemcpy(ZR,(double *) d_VR, sizeof(double)*nr, cudaMemcpyDeviceToHost));
    fwrite(Time,sizeof(double),1,avgs);
    fwrite(ZR,sizeof(double),nr,avgs);
    
    subtractReal<<<nblocks,nthreads>>>(d_UR,d_Urms);
    
    (cudaMemcpy(ZR,(double *) d_Vrms, sizeof(double)*nr, cudaMemcpyDeviceToHost));
    fwrite(Time,sizeof(double),1,avgs);
    fwrite(ZR,sizeof(double),nr,avgs);
    
    (cudaMemcpy(ZR,(double *) d_UR, sizeof(double)*nr, cudaMemcpyDeviceToHost));
    fwrite(Time,sizeof(double),1,avgs);
    fwrite(ZR,sizeof(double),nr,avgs);
    
    fflush(avgs);
    
    fclose(avgs);
    fclose(stats);
    fclose(points);

    (cudaFree(d_Urms));
    (cudaFree(d_Vrms));
    //    free(Zstore);
    free(ZKtavg);
  }


  fclose(vort);
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
  (cufftDestroy(PlanZ2D));
  (cufftDestroy(PlanBatchD2Z));


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
