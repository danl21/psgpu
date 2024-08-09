#define TWPI 6.283185307179586
#define MY 0
////////////////////////////////////////////////////////////////////////
// Define constant device variables
////////////////////////////////////////////////////////////////////////
__constant__ int d_IKTX, d_IKTY, d_KTY, d_NX, d_NY,d_NX2,d_KFX,d_KFY,d_OR,d_OK,d_O2,d_nsx,d_nsy;
__constant__ double d_IN,d_DELT,d_NU;

// Define constant host variables
////////////////////////////////////////////////////////////////////////

// Threads should be a multiple of warp size (32) and maximise warps per multiproc to hide register latency (>192)
// Blocks should be a multiple of multiprocs (14) and >100 for scalability
const int nthreads=512;
const int nblocks=140;

int nkt,iktx,ikty,kty,nr,nx,ny,nx2,n2,nOut,nStore,izkout,RecOUTflag,minSTOREouter,iStart,tavgCount;
int np[2];
double v2,in,freq2,tme,tstart,delt,ResidualThreshold,energy,diss,input,ampfor,section,ZPP1,ZPP2;

const int storeSize = 500;
const int normWidth= 42;
const int nWidth2 = normWidth;
const int nNorm=2*normWidth*nWidth2;
const double twopi = 4.0*asin(1.0);
const int NSHIFTX =30;
const int NSHIFTY =4;

// recurrence variables
int iNorm[nNorm];
double shiftX[NSHIFTX],shiftY[NSHIFTY],normZStore[storeSize];
double minRESIDUALouter, minTIMEouter,RecPERIOD,minTIMEinner,normZ;
double minSHIFTXouter,minSHIFTYouter,minSHIFTXinner,minSHIFTYinner;
double *d_kx,*d_ky;
FILE * UPOfile;
FILE * UPOZK;
FILE * stats;
FILE * points;
FILE * avgs;
FILE * vort;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//******************  GPU KERNELS    ***************************
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void setFFN(cufftDoubleComplex *V, cufftDoubleComplex *F1, int *ikF){
  int ik = threadIdx.x + blockIdx.x*blockDim.x;

  while(ik < d_O2){
    // Map -Ky values to -Ky+NY and fill trucated wave numbers with zeros
    // Note ikF is set with indexing for this on the fortran side
      int ikk = ikF[ik];
      if(ikk < 0){
	F1[ik].x = 0.0;
	F1[ik].y = 0.0;
      }else{
	F1[ik].x = V[ikk].x;
	F1[ik].y = V[ikk].y;
      }
    
    ik+= blockDim.x*gridDim.x;
  }
}


__global__ void avgVelocity(cufftDoubleReal *UR,cufftDoubleReal *VR,cufftDoubleReal *Urms,cufftDoubleReal *Vrms, int *LL, int tavgCount){
  // Kernel to set the spectral velocity coefficients UK and VK for first step
  int ik = threadIdx.x + blockIdx.x*blockDim.x;

  while(ik < d_OR){

    Urms[ik] = (Urms[ik]*tavgCount + UR[ik]*UR[ik])/(tavgCount+1.0);
    Vrms[ik] = (Vrms[ik]*tavgCount + VR[ik]*VR[ik])/(tavgCount+1.0);

    ik += blockDim.x*gridDim.x;
  }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\
///////////

__global__ void conj(cufftDoubleComplex *ZK, double *kx, double *ky){
  // Kernel to set the spectral velocity coefficients UK and VK for first step
  int ik = threadIdx.x + blockIdx.x*blockDim.x;

  while(ik < d_OK){
    if( kx[ik] == 0.0 && ky[ik] < 0.0){
      int iky = ik/d_IKTX;
      int ikk = (d_IKTY-iky-1)*d_IKTX;

      ZK[ik].x =  ZK[ikk].x;
      ZK[ik].y = -ZK[ikk].y;
    }
    ik += blockDim.x*gridDim.x;
  }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\
///////////

__global__ void subLam(cufftDoubleComplex *ZK,double v2){
  // Kernel to set the spectral velocity coefficients UK and VK for first step
  int ik = (d_KTY+4)*d_IKTX;
  ZK[ik].x += 0.125/v2;

}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\
/////////

__global__ void setVelocity(cufftDoubleComplex *ZK,cufftDoubleComplex *Z0,cufftDoubleComplex *UK,cufftDoubleComplex *VK,double *kx,double *ky, int *LL){
  // Kernel to set the spectral velocity coefficients UK and VK for first step
  int ik = threadIdx.x + blockIdx.x*blockDim.x;

  while(ik < d_OK){

  // Apply mask
    if(LL[ik]==1){
      double kkx = kx[ik];
      double kky = ky[ik];
      
      double wk =kkx*kkx + kky*kky;
      double wki = 1.0/(max(wk,0.001));

      UK[ik].x = -kky*ZK[ik].y*wki;
      UK[ik].y =  kky*ZK[ik].x*wki;
      
      VK[ik].x =  kkx*ZK[ik].y*wki;
      VK[ik].y = -kkx*ZK[ik].x*wki;

      Z0[ik].x = ZK[ik].x;
      Z0[ik].y = ZK[ik].y;
      
    }else{
      UK[ik].x = 0.0;
      UK[ik].y = 0.0;
      
      VK[ik].x = 0.0;
      VK[ik].y = 0.0;

      Z0[ik].x = 0;
      Z0[ik].y = 0;
    }
    ik += blockDim.x*gridDim.x;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void normFF1(cufftDoubleComplex *F2, cufftDoubleComplex *F1, int *ikN){
  int ik = threadIdx.x + blockIdx.x*blockDim.x;
  
  //normalise D2Z FFT output
  while(ik < 2*d_OK){

    if(ik < d_OK){
      int ikk = ikN[ik];
      F2[ik].x = F1[ikk].x*d_IN;
      F2[ik].y = F1[ikk].y*d_IN;
    }else{
      int ikk = ikN[ik-d_OK]+d_O2;
      F2[ik].x = F1[ikk].x*d_IN;
      F2[ik].y = F1[ikk].y*d_IN;
    }
    ik += blockDim.x*gridDim.x;
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void normFFZ(cufftDoubleComplex *F2, cufftDoubleComplex *F1, int *ikN){
  int ik = threadIdx.x + blockIdx.x*blockDim.x;

  //normalise D2Z FFT output
  while(ik < d_OK){

    int ikk = ikN[ik];
    F2[ik].x = F1[ikk].x*d_IN;
    F2[ik].y = F1[ikk].y*d_IN;
    ik += blockDim.x*gridDim.x;
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void multReal(cufftDoubleReal* Z,cufftDoubleReal* U,cufftDoubleReal* V, cufftDoubleReal* NZ){
  int ik = threadIdx.x + blockIdx.x*blockDim.x;
  // do physical space multiplication (cross product U X Z)

  while(ik < 2*d_OR){
    if(ik < d_OR){
      NZ[ik] = V[ik]*Z[ik];
    }else{
      NZ[ik] = U[ik-d_OR]*Z[ik-d_OR];
    }
    ik += blockDim.x*gridDim.x;
  }

  
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void zeroComplex(cufftDoubleComplex *F1){
  int ik = threadIdx.x + blockIdx.x*blockDim.x;

  while(ik < d_OK){
    F1[ik].x = 0.0;
    F1[ik].y = 0.0;

    F1[ik+d_OK].x = 0.0;
    F1[ik+d_OK].y = 0.0;

    ik += blockDim.x*gridDim.x;
  }
  
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void zeroReal(cufftDoubleReal *F1){
  int ik = threadIdx.x + blockIdx.x*blockDim.x;

  while(ik < d_OR){
    F1[ik] = 0.0;

    ik += blockDim.x*gridDim.x;
  }
  
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void subtractReal(cufftDoubleReal *F1,cufftDoubleReal *F2){
  int ik = threadIdx.x + blockIdx.x*blockDim.x;

  while(ik < d_OR){
    F1[ik] = F2[ik]-F1[ik]*F1[ik];

    ik += blockDim.x*gridDim.x;
  }
  
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void Step1(cufftDoubleComplex *Z,cufftDoubleComplex *Z0,cufftDoubleComplex *Z1,cufftDoubleComplex *UK,cufftDoubleComplex *VK,cufftDoubleComplex *NZK,double *kx,double *ky,int *LL){
  int ik = threadIdx.x + blockIdx.x*blockDim.x;
  
  while(ik < d_OK){
    if(LL[ik]==1){

      double ZX,ZY;
      double kky = ky[ik];
      double kkx = kx[ik];
      double wk = kkx*kkx + kky*kky;
      double wki = 1.0/(max(wk,0.001));

      // Do final part of 'convol' (CURL)
      double NZKX =  kkx*NZK[ik+d_OK].y + kky*NZK[ik].y;
      double NZKY = -kkx*NZK[ik+d_OK].x - kky*NZK[ik].x;

      //check if this is a forcing wave number and make appropriate adjustments
      if(kkx ==0.0 && abs(kky) == 4.0){
	NZKX -= 2.0;
      }

      ZX = d_DELT*NZKX;
      ZY = d_DELT*NZKY;
      
      // Update vorticity array
      Z1[ik].x = ZX;
      Z1[ik].y = ZY;

      Z0[ik].x = Z[ik].x + 0.5*ZX;
      Z0[ik].y = Z[ik].y + 0.5*ZY;

      // Update spectral velocity coeffs
      UK[ik].x = -kky*Z0[ik].y*wki;
      UK[ik].y =  kky*Z0[ik].x*wki;

      VK[ik].x =  kkx*Z0[ik].y*wki;
      VK[ik].y = -kkx*Z0[ik].x*wki;
      
    }else{
      UK[ik].x = 0.0;
      UK[ik].y = 0.0;

      VK[ik].x = 0.0;
      VK[ik].y = 0.0;

      Z[ik].x = 0.0;
      Z[ik].y = 0.0;
      
    }

    // reset this array
    NZK[ik].x = 0.0;
    NZK[ik].y = 0.0;

    NZK[ik+d_OK].x = 0.0;
    NZK[ik+d_OK].y = 0.0;

    ik+=blockDim.x*gridDim.x;
    
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void Step2(cufftDoubleComplex *Z,cufftDoubleComplex *Z0,cufftDoubleComplex *Z1,cufftDoubleComplex *UK,cufftDoubleComplex *VK,cufftDoubleComplex *NZK,double *kx,double *ky,int *LL){
  int ik = threadIdx.x + blockIdx.x*blockDim.x;
  
  while(ik < d_OK){
    if(LL[ik]==1){

      double ZX,ZY;
      double kky = ky[ik];
      double kkx = kx[ik];
      double wk =kkx*kkx + kky*kky;
      double wki = 1.0/(max(wk,0.001));

      // Do final part of 'convol' (CURL)
      double NZKX =  kkx*NZK[ik+d_OK].y + kky*NZK[ik].y;
      double NZKY = -kkx*NZK[ik+d_OK].x - kky*NZK[ik].x;

      //check if this is a forcing wave number and make appropriate adjustments
      if(kkx ==0.0 && abs(kky) == 4.0){
	NZKX -= 2.0;
      }

      ZX = d_DELT*NZKX;
      ZY = d_DELT*NZKY;
      
      // Update vorticity array
      Z1[ik].x = ZX;
      Z1[ik].y = ZY;

      Z0[ik].x = Z[ik].x + ZX;
      Z0[ik].y = Z[ik].y + ZY;

      // Update spectral velocity coeffs
      UK[ik].x = -kky*Z0[ik].y*wki;
      UK[ik].y =  kky*Z0[ik].x*wki;

      VK[ik].x =  kkx*Z0[ik].y*wki;
      VK[ik].y = -kkx*Z0[ik].x*wki;
      
    }else{
      UK[ik].x = 0.0;
      UK[ik].y = 0.0;

      VK[ik].x = 0.0;
      VK[ik].y = 0.0;

      Z[ik].x = 0.0;
      Z[ik].y = 0.0;
      
    }

    // reset this array
    NZK[ik].x = 0.0;
    NZK[ik].y = 0.0;

    NZK[ik+d_OK].x = 0.0;
    NZK[ik+d_OK].y = 0.0;

    ik+=blockDim.x*gridDim.x;
    
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void Step(cufftDoubleComplex *Z,cufftDoubleComplex *Z0,cufftDoubleComplex *Z1,cufftDoubleComplex *Z2,cufftDoubleComplex *Z3,cufftDoubleComplex *UK,cufftDoubleComplex *VK,cufftDoubleComplex *NZK,double *kx,double *ky,int *LL){
  int ik = threadIdx.x + blockIdx.x*blockDim.x;
  
  while(ik < d_OK){
    if(LL[ik]==1){

      double ZX,ZY;
      double kky = ky[ik];
      double kkx = kx[ik];
      double wk =kkx*kkx + kky*kky;
      double wki = 1.0/(max(wk,0.001));

      // Do final part of 'convol' (CURL)
      double NZKX =  kkx*NZK[ik+d_OK].y + kky*NZK[ik].y;
      double NZKY = -kkx*NZK[ik+d_OK].x - kky*NZK[ik].x;

      double damp1 = 1-0.5*d_DELT*d_NU*wk;
      double damp2 = 1./(1+0.5*d_DELT*d_NU*wk);
      //check if this is a forcing wave number and make appropriate adjustments
      if(kkx ==0.0 && abs(kky) == 4.0){
	NZKX -= 2.0;
      }

      ZX = d_DELT*NZKX;
      ZY = d_DELT*NZKY;
      
      // Update vorticity array
      Z[ik].x = damp2*(damp1*Z[ik].x + (Z1[ik].x+2*Z2[ik].x+2*Z3[ik].x+ZX)/6.);
      Z[ik].y = damp2*(damp1*Z[ik].y + (Z1[ik].y+2*Z2[ik].y+2*Z3[ik].y+ZY)/6.);

      Z0[ik].x = Z[ik].x;
      Z0[ik].y = Z[ik].y;

      // Update spectral velocity coeffs
      UK[ik].x = -kky*Z0[ik].y*wki;
      UK[ik].y =  kky*Z0[ik].x*wki;

      VK[ik].x =  kkx*Z0[ik].y*wki;
      VK[ik].y = -kkx*Z0[ik].x*wki;
      
    }else{
      UK[ik].x = 0.0;
      UK[ik].y = 0.0;

      VK[ik].x = 0.0;
      VK[ik].y = 0.0;

      Z[ik].x = 0.0;
      Z[ik].y = 0.0;
      
    }

    // reset this array
    NZK[ik].x = 0.0;
    NZK[ik].y = 0.0;

    NZK[ik+d_OK].x = 0.0;
    NZK[ik+d_OK].y = 0.0;

    ik+=blockDim.x*gridDim.x;
    
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__host__ __device__ double carg(const cuDoubleComplex& z) {return atan2(cuCimag(z),   cuCreal(z));}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void mode_shift(cufftDoubleComplex *Z1,cufftDoubleComplex *Z2,double *SS,double *S1,double *kx, int *LL){
  // Kernel to set the spectral velocity coefficients UK and VK for first step                                       
  int ik = threadIdx.x + blockIdx.x*blockDim.x;
  cuDoubleComplex z1,z2,epsC = make_cuDoubleComplex(1.e-32,1.e-32);
  double eps  = 1.e-32;

  while(ik < d_OK){

  // Apply mask                                                                                                      
    if(LL[ik]==1 && kx[ik]==1){

      z1 = Z1[ik];
      z2 = Z2[ik];
      S1[ik] = carg(cuCdiv(z1,cuCadd(z2,epsC)))/(kx[ik]+eps);
      SS[ik] = S1[ik];

    }else if(LL[ik]==1 && kx[ik]!=0){
      z1 = Z1[ik];
      z2 = Z2[ik];
      SS[ik] = carg(cuCdiv(z1,cuCadd(z2,epsC)))/(kx[ik]+eps);

      S1[ik] = 0.0;
    }else{
      SS[ik] = 0.0;
      S1[ik] = 0.0;
    }
    ik += blockDim.x*gridDim.x;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\
///////////////                                                                                                      

__global__ void shift_branch(double *SS,double *S1, double *Smean, double *kx,int *LL){
  // Kernel to set the spectral velocity coefficients UK and VK for first step                                       
  int ik = threadIdx.x + blockIdx.x*blockDim.x;
  double eps  = 1.e-32;

  while(ik < d_OK){

  // Apply mask                                                                                                      
    if(LL[ik]==1 && kx[ik]!=0){
      double temp = (*Smean-SS[ik])*kx[ik];
      int nN = round(temp/TWPI);
      S1[ik] = (SS[ik]*kx[ik]+nN*TWPI)/(kx[ik]+eps);
    }else{
      SS[ik] = 0.0;
      S1[ik] = 0.0;
    }
    ik += blockDim.x*gridDim.x;
  }
}

__global__ void symmZ(cufftDoubleComplex *Z,cufftDoubleComplex *Z0,double *kx,double *ky, int *LL){
  // Kernel to set the spectral velocity coefficients UK and VK for first step                                                                                               
  int ik = threadIdx.x + blockIdx.x*blockDim.x;

  while(ik < d_OK){

  // Apply mask                                                                                                                                                              
    if(LL[ik]==1){
      double kkx = kx[ik];
      double kky = ky[ik];

      double wk =kkx*kkx + kky*kky;
      double wki = 1.0/(max(wk,0.001));

      int iky = ik/d_IKTX-d_KTY;
      int ikx = ik - (iky+d_KTY)*d_IKTX;
      int iky2 = -iky;
      int ikk = d_IKTX*(iky2+d_KTY) + ikx;
      double ss,cs,shift=-kkx*(0.5*TWPI)-kky*MY*TWPI*0.5;
      sincos(shift,&ss,&cs);

      if(kx[ik]==0.0){
        Z0[ik].x = -(cs*Z[ik].x+ss*Z[ik].y);
        Z0[ik].y =  (cs*Z[ik].y-ss*Z[ik].x);
      }else{
        Z0[ikk].x = -(cs*Z[ik].x+ss*Z[ik].y);
        Z0[ikk].y = -(cs*Z[ik].y-ss*Z[ik].x);
      }
    }else{
      Z0[ik].x = 0;
      Z0[ik].y = 0;
    }
    ik += blockDim.x*gridDim.x;
  }
}
