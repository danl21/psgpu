///////////////////////////////////////////////////////////////////////////////

extern "C" void set_gpu_(int *go)
{
  int            dev, deviceCount;
  cudaDeviceProp devProp;

  (cudaGetDeviceCount(&deviceCount));
  if (deviceCount == 0) {
    printf("cutil error: no devices supporting CUDA\n");
    exit(-1);
  }

  (cudaGetDevice(&dev));
  //  (cudaSetDevice(*go));
  (cudaGetDeviceProperties(&devProp,dev));
  printf("\n Using CUDA device %d: %s\n\n", dev,devProp.name);


}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\
////////////////////////

extern "C" void KR_FFT_ALL(cufftDoubleComplex *d_xii,
			   cufftDoubleComplex *d_eta,
			   cufftDoubleComplex *d_zet, 
			   cufftDoubleComplex *d_UK,
			   cufftDoubleComplex *d_VK,
			   cufftDoubleComplex *d_WK, 
			   cufftDoubleReal*d_xir, 
			   cufftDoubleReal*d_etr, 
			   cufftDoubleReal*d_ztr, 
			   cufftDoubleReal*d_UR, 
			   cufftDoubleReal*d_VR,
			   cufftDoubleReal*d_WR,
			   cufftHandle PlanZ2D,
			   int *d_ikF, 
			   int n2){

  cufftDoubleComplex *d_FF,*d_tmp;

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
  (cufftExecZ2D(PlanZ2D,d_FF,d_WR));  // Do V FFT

  (cudaFree(d_FF));
  (cudaFree(d_tmp));
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\
////////////////////////

extern "C" void RK_FFT(cufftDoubleComplex *d_NXK,
		       cufftDoubleComplex *d_NYK,
		       cufftDoubleComplex *d_NZK, 
		       cufftDoubleReal*d_XR,
		       cufftDoubleReal*d_YR,
		       cufftDoubleReal*d_ZR,
		       cufftHandle PlanD2Z,
		       int *d_ikN,
		       int n2){
  
  cufftDoubleComplex *d_FF;

  //Doing a batch of 2 2D FFTs so that now NZK stores the spectral UK and NZK from the old convol in that order.

  (cudaMalloc((void**)&d_FF,sizeof(cufftDoubleComplex)*(n2)));

  (cufftExecD2Z(PlanD2Z,d_XR,d_FF));

  normFF1<<<nblocks,nthreads>>>(d_NXK,d_FF,d_ikN);  // normalise output

  (cufftExecD2Z(PlanD2Z,d_YR,d_FF));

  normFF1<<<nblocks,nthreads>>>(d_NYK,d_FF,d_ikN);  // normalise output

  (cufftExecD2Z(PlanD2Z,d_ZR,d_FF));

  normFF1<<<nblocks,nthreads>>>(d_NZK,d_FF,d_ikN);  // normalise output

  (cudaFree(d_FF));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

extern "C" void computeStats(double2 *XII,double2 *ETA,double2 *ZET,double2 *XIItavg,double2 *ETAtavg,double2 *ZETtavg,double *kx, double *ky,double *kz,int *LL){

  double eturb = 0.0;
  double enstrophy = 0.0;
  double XIa = 0;
  double ETa = 0;
  double ZTa = 0;
  double Ua = 0;
  double Va = 0;
  double Wa = 0;

  energy = 0.0;
  for(int i=0; i<nkt; i++){
    if(LL[i] == 1){
    double wk = max(kx[i]*kx[i] + ky[i]*ky[i]+kz[i]*kz[i],0.001);

    double VZX = XII[i].x*XII[i].x + XII[i].y*XII[i].y;
    double VZY = ETA[i].x*ETA[i].x + ETA[i].y*ETA[i].y;
    double VZZ = ZET[i].x*ZET[i].x + ZET[i].y*ZET[i].y;
    double VZ = VZX+VZY+VZZ;


    energy += VZ/wk;
    enstrophy += VZ;
    
    XIa+=VZX;
    ETa+=VZY;
    ZTa+=VZZ;
    Ua += pow((ky[i]*ZET[i].y-kz[i]*ETA[i].y)/wk,2) + pow((kz[i]*ETA[i].x-ky[i]*ZET[i].x)/wk,2);
    Va += pow((kz[i]*XII[i].y-kx[i]*ZET[i].y)/wk,2) + pow((kx[i]*ZET[i].x-kz[i]*XII[i].x)/wk,2);
    Wa += pow((kx[i]*ETA[i].y-ky[i]*XII[i].y)/wk,2) + pow((ky[i]*XII[i].x-kx[i]*ETA[i].x)/wk,2);

    XIItavg[i].x = (XIItavg[i].x*tavgCount + XII[i].x)/(tavgCount+1.0);
    XIItavg[i].y = (XIItavg[i].y*tavgCount + XII[i].y)/(tavgCount+1.0);
    ETAtavg[i].x = (ETAtavg[i].x*tavgCount + ETA[i].x)/(tavgCount+1.0);
    ETAtavg[i].y = (ETAtavg[i].y*tavgCount + ETA[i].y)/(tavgCount+1.0);
    ZETtavg[i].x = (ZETtavg[i].x*tavgCount + ZET[i].x)/(tavgCount+1.0);
    ZETtavg[i].y = (ZETtavg[i].y*tavgCount + ZET[i].y)/(tavgCount+1.0);

    VZX = XIItavg[i].x*XIItavg[i].x + XIItavg[i].y*XIItavg[i].y;
    VZY = ETAtavg[i].x*ETAtavg[i].x + ETAtavg[i].y*ETAtavg[i].y;
    VZZ = ZETtavg[i].x*ZETtavg[i].x + ZETtavg[i].y*ZETtavg[i].y;
    VZ = VZX+VZY+VZZ;

    eturb += VZ/wk;
    }
  }
  enerMax = max(energy,enerMax);
  enerMin = min(energy,enerMin);
  enerAvg = (enerAvg*tavgCount + energy)/(tavgCount+1.0);

  int istat= int(NT/nOut)-1;
  Energy[istat]= energy;

  tavgCount++;
  fprintf(stats,"%e %e %e \n",tme,energy,eturb,enstrophy);
  fflush(stats);
  fprintf(engy,"%e %e %e %e %e %e %e \n",tme,Ua,Va,Wa,XIa,ETa,ZTa);
  fflush(engy);
  
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

