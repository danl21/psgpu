////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// THIS IS A SOURCE FILE FOR THE 3D_PSPGU CODE 
// IT CONTAINS THE DEVICE KERNELS USED FOR
// TIMESTEPPING THE PROGNOSTIC EQS AND SETTING 
// SOME GLOBAL VARIABLES.
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// Some compiler macro constants (host and/or device)
#define sixth (1.0/6.0)
#define third (1.0/3.0)
#define storeSize 500
#define normWidth 8
#define NSHIFTX 64
#define NSHIFTZ 1
#define NSHIFTY 1
////////////////////////////////////////////////////////////////////////
// Define constant device variables
////////////////////////////////////////////////////////////////////////
__constant__ int d_IKTX,d_IKTY,d_IKTZ,d_KTZ,d_NX, d_NY,d_NZ,d_NX2,d_OR,d_OK,d_O2;
__constant__ double d_IN,d_AMPFOR,d_DELT,d_v2,d_KFY,d_Ri,d_Sc,d_SinTh,d_CosTh;

// Define constant host variables
////////////////////////////////////////////////////////////////////////

// Threads should be a multiple of warp size (32) and maximise warps per multiproc to hide register latency (>192)
// Blocks should be a multiple of multiprocs (14) and >100 for scalability
const int nthreads=512;
const int nblocks=140;
const double twopi = 4.0*asin(1.0);
const double small = 1e-10;
const int nNorm = 4*normWidth*normWidth*normWidth;
////////////////////////////////////////////////////////////////////////
// host variables
int nkt,iktx,ikty,iktz,nr,nx,ny,nz,nx2,n2,nOut,nStore,izkout,iStart,tavgCount,NT,rank,countStats;
double v2,in,freq2,tme,tstart,delt,ResidualThreshold,energy,diss,input,ampfor,kfy,pe,ped,bflux,sinTh,cosTh,Umax,ri,dtarget;
double peAvg,enerAvg, dissAvg, inputAvg,eturbAvg,dturbAvg,mixAvg,pedAvg,bfluxAvg,enerMin,dissMin,dissMax,enerMax,inputMin,inputMax;

// recurrence variables
int iNorm[nNorm],RecOUTflag,minSTOREouter;
double shiftX[NSHIFTX],shiftY[NSHIFTY],shiftZ[NSHIFTZ],normZStore[storeSize];
double minRESIDUALouter, minTIMEouter,RecPERIOD,normZ;
double minSHIFTXouter,minSHIFTYouter,minSHIFTZouter,minSHIFTXinner,minSHIFTYinner,minSHIFTZinner;
// device pointers for wavevector
double *d_kx,*d_ky,*d_kz;
// host pointers for statistics
double *Energy,*Diss,*Input;
// output streams
FILE *UPOfile,*UPOZK,*sweep,*ZK,*stats,*engy,*avgs,*means,*pdfs,*spec,*RiG,*bfile,*xiiStream,*etaStream,*zetStream,*rhoStream;

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
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\

__global__ void setFFK(cufftDoubleComplex *V, cufftDoubleComplex *F1,double *kk, int *ikF){
  int ik = threadIdx.x + blockIdx.x*blockDim.x;

  while(ik < d_O2){
    // Map -Ky values to -Ky+NY and fill trucated wave numbers with zeros
    // Note ikF is set with indexing for this on the fortran side
      int ikk = ikF[ik];
      if(ikk < 0){
      	F1[ik].x = 0.0;
      	F1[ik].y = 0.0;
      }else{
	// compute a gradient at the same time...
	F1[ik].x = -kk[ikk]*V[ikk].y;
	F1[ik].y =  kk[ikk]*V[ikk].x;
      }
    
    ik+= blockDim.x*gridDim.x;
  }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\
//Kernel to average the rms velocity
__global__ void avgVelocity(cufftDoubleReal *UR,cufftDoubleReal *VR,cufftDoubleReal *WR,cufftDoubleReal *Urms,cufftDoubleReal *Vrms,cufftDoubleReal *Wrms, int tavgCount){
  int ik = threadIdx.x + blockIdx.x*blockDim.x;

  while(ik < d_OR){

    Urms[ik] = (Urms[ik]*tavgCount + UR[ik]*UR[ik])/(tavgCount+1.0);
    Vrms[ik] = (Vrms[ik]*tavgCount + VR[ik]*VR[ik])/(tavgCount+1.0);
    Wrms[ik] = (Wrms[ik]*tavgCount + WR[ik]*WR[ik])/(tavgCount+1.0);

    ik += blockDim.x*gridDim.x;
  }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\
//////////
// belt and braces check that arrays are full of 0.0
__global__ void zeroReal(cufftDoubleReal *F1){
  int ik = threadIdx.x + blockIdx.x*blockDim.x;

  while(ik < d_OR){
    F1[ik] = 0.0;

    ik += blockDim.x*gridDim.x;
  }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\
//////////
// kernel to subtract two real space arrays
__global__ void subtractReal(cufftDoubleReal *F1,cufftDoubleReal *F2){
  int ik = threadIdx.x + blockIdx.x*blockDim.x;

  while(ik < d_OR){
    F1[ik] = F2[ik]-F1[ik]*F1[ik];

    ik += blockDim.x*gridDim.x;
  }

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void conj(cufftDoubleComplex *ZK,cufftDoubleComplex *ZZ, double *kx){
  int ik = threadIdx.x + blockIdx.x*blockDim.x;
  // Enforce conjugate symmetry plane (this may be redundant for some versions of CUFFT)
  while(ik < d_OK){ 
    if( kx[ik] == 0.0 ){
      int ikt = int(ik/d_IKTX);
      int ikz = int(floor(double(ikt / d_IKTY)));
      int iky = int(ikt - ikz*d_IKTY);

      if( iky == d_KTZ ){
	if( ikz == d_KTZ){
	  ZZ[ik].x = 0.0;
	  ZZ[ik].y = 0.0;
	}else if(ikz < d_KTZ){
	int ikzz = d_IKTZ-1 - ikz;

	int ikk = (ikzz*d_IKTY + iky)*d_IKTX;
	double ZKX = ZK[ikk].x; 
	double ZKY =-ZK[ikk].y;
	
	ZZ[ik].x =  ZKX;
	ZZ[ik].y =  ZKY;	  
	}else{
	ZZ[ik].x = ZK[ik].x;
	ZZ[ik].y = ZK[ik].y;
	}
      }else if( iky < d_KTZ){
	int ikky = d_IKTY-1 -iky;
	
	ikz = d_IKTZ-1 - ikz;

	int ikk = (ikz*d_IKTY + ikky)*d_IKTX;
	
	double ZKX = ZK[ikk].x; 
	double ZKY =-ZK[ikk].y;
	
	ZZ[ik].x =  ZKX;
	ZZ[ik].y =  ZKY;
	}else{

	ZZ[ik].x = ZK[ik].x;
	ZZ[ik].y = ZK[ik].y;
      }
      }else{
	ZZ[ik].x = ZK[ik].x;
	ZZ[ik].y = ZK[ik].y;
    }
      ik += blockDim.x*gridDim.x;
  }

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void setVelocity(cufftDoubleComplex *xii,cufftDoubleComplex *eta,cufftDoubleComplex *zet,cufftDoubleComplex *rho,cufftDoubleComplex *x0,cufftDoubleComplex *e0,cufftDoubleComplex *z0,cufftDoubleComplex *r0,cufftDoubleComplex *UK,cufftDoubleComplex *VK,cufftDoubleComplex *WK,double *kx,double *ky,double *kz,int *LL){
  // Kernel to set the spectral velocity coefficients UK and VK for first step
  int ik = threadIdx.x + blockIdx.x*blockDim.x;

  while(ik < d_OK){
      double kkx = kx[ik];
      double kky = ky[ik];
      double kkz = kz[ik];

      double xiX = xii[ik].x;
      double xiY = xii[ik].y;
      
      double etX = eta[ik].x;
      double etY = eta[ik].y;
      
      double ztX = zet[ik].x;
      double ztY = zet[ik].y;
      

      double wk =kkx*kkx + kky*kky +kkz*kkz;
      double wki = 1.0/(max(wk,0.001));

	// set mask this way to minimise warp divergence
      if(LL[ik] == 0) wki = 0.0;

      // Update spectral velocity coeffs
      UK[ik].x = -(kky*ztY-kkz*etY)*wki;
      UK[ik].y = -(kkz*etX-kky*ztX)*wki;

      VK[ik].x = -(kkz*xiY-kkx*ztY)*wki;
      VK[ik].y = -(kkx*ztX-kkz*xiX)*wki;

      WK[ik].x = -(kkx*etY-kky*xiY)*wki;
      WK[ik].y = -(kky*xiX-kkx*etX)*wki;

      x0[ik].x = xiX;
      x0[ik].y = xiY;
      e0[ik].x = etX;
      e0[ik].y = etY;
      z0[ik].x = ztX;
      z0[ik].y = ztY;
      r0[ik].x = rho[ik].x;
      r0[ik].y = rho[ik].y;

      ik += blockDim.x*gridDim.x;
  }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void zeroVorticity(cufftDoubleComplex *xii,cufftDoubleComplex *eta,cufftDoubleComplex *zet){
  // Kernel to set the spectral vorticity to zero initially
  int ik = threadIdx.x + blockIdx.x*blockDim.x;

  while(ik < d_OK){

      xii[ik].x = 0.0;
      xii[ik].y = 0.0;
      
      eta[ik].x = 0.0;
      eta[ik].y = 0.0;

      zet[ik].x = 0.0;
      zet[ik].y = 0.0;


      ik += blockDim.x*gridDim.x;
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void setVorticity(cufftDoubleComplex *xii,cufftDoubleComplex *eta,cufftDoubleComplex *zet,double *kx,double *ky,double *kz){
  // Kernel to set the spectral vorticity for initial condition
  int ik = threadIdx.x + blockIdx.x*blockDim.x;

  while(ik < d_OK){

  if(kx[ik]==0.0 && ky[ik]== 1.0 && kz[ik]== 1.0) xii[ik].x -= 0.5;
  if(kx[ik]==0.0 && ky[ik]== -1.0 && kz[ik]==-1.0) xii[ik].x -= 0.5;
  if(kx[ik]==0.0 && ky[ik]== 1.0 && kz[ik]==-1.0) xii[ik].x +=  0.5;
  if(kx[ik]==0.0 && ky[ik]== -1.0 && kz[ik]==1.0) xii[ik].x +=  0.5;

  ik += blockDim.x*gridDim.x;
  }

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void normFF1(cufftDoubleComplex *F2, cufftDoubleComplex *F1, int *ikN,int *LL){
  int ik = threadIdx.x + blockIdx.x*blockDim.x;
  
  //normalise D2Z FFT output d_IN holds the array size
  while(ik < d_OK){

    int ikk = ikN[ik];
    F2[ik].x = F1[ikk].x*d_IN;
    F2[ik].y = F1[ikk].y*d_IN;
    ik += blockDim.x*gridDim.x;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void multReal(cufftDoubleReal* xii,cufftDoubleReal* eta,cufftDoubleReal* zet,
			 cufftDoubleReal* U,cufftDoubleReal* V,cufftDoubleReal* W,
			 cufftDoubleReal* rx,cufftDoubleReal* ry, cufftDoubleReal* rz){
  int ik = threadIdx.x + blockIdx.x*blockDim.x;
  // do physical space multiplication (cross product U X Z)
  // Note careful use of local variables to minimise cache use

  while(ik < d_OR){
    double TU,TV,TW,TE,TZ,TX,NX,NY;
    TV = V[ik];
    TU = U[ik];
    TW = W[ik];

    TX = xii[ik];
    TE = eta[ik];
    TZ = zet[ik];

    NX = TV*TZ-TW*TE;
    NY = TW*TX-TU*TZ;
    TZ = TU*TE-TV*TX;
    
    U[ik]=-NX;
    V[ik]=-NY;
    W[ik]=-TZ;
    rx[ik]= TU*rx[ik]+TV*ry[ik]+TW*rz[ik];

    ik += blockDim.x*gridDim.x;
  }

  
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void Step1(cufftDoubleComplex *xii,cufftDoubleComplex *eta,cufftDoubleComplex *zet,cufftDoubleComplex *rho,
		      cufftDoubleComplex *x0,cufftDoubleComplex *e0,cufftDoubleComplex *z0,cufftDoubleComplex *r0,
		      cufftDoubleComplex *x1,cufftDoubleComplex *e1,cufftDoubleComplex *z1,cufftDoubleComplex *r1,
		      cufftDoubleComplex *UK,cufftDoubleComplex *VK,cufftDoubleComplex *WK,cufftDoubleComplex *RK,
		      double *kx,double *ky,double *kz, int *LL){
  int ik = threadIdx.x + blockIdx.x*blockDim.x;
  // This kernel does the first type of sub-step for RK4
  // Most quantities are loaded into local variables to help caching 
  while(ik < d_OK){
      double kky = ky[ik];
      double kkx = kx[ik];
      double kkz = kz[ik];
      double UKX = UK[ik].x;
      double UKY = UK[ik].y;
      double VKX = VK[ik].x;
      double VKY = VK[ik].y;
      double WKX = WK[ik].x;
      double WKY = WK[ik].y;
      double ROX = r0[ik].x;
      double ROY = r0[ik].y;

      double wk  = kkx*kkx + kky*kky +kkz*kkz;
      double wki = 1.0/(max(wk,0.001));

      // gradients of density
      double rXR = -kkx*ROY;
      double rXC =  kkx*ROX;
      double rYR = -kky*ROY;
      double rYC =  kky*ROX;
      double rZR = -kkz*ROY;
      double rZC =  kkz*ROX;
 
      // dealiasing mask (2/3rds or Hou filter, this choice is important for GMRES too!)
      double LC = LL[ik];
      //double LC = exp(-36.0*(pow(2*sqrt(wk)/double(d_NX),36.0)));

      // Do final part of 'convol' (CURL) and add buoyancy term projected due to angle of g
      double xiX =  LC*(kky*WKY - kkz*VKY- d_Ri*(rYR*d_SinTh - rZR*d_CosTh));
      double xiY =  LC*(kkz*VKX - kky*WKX- d_Ri*(rYC*d_SinTh - rZC*d_CosTh));

      double etX =  LC*(kkz*UKY - kkx*WKY+ d_Ri*rXR*d_SinTh);
      double etY =  LC*(kkx*WKX - kkz*UKX+ d_Ri*rXC*d_SinTh);

      double ztX =  LC*(kkx*VKY - kky*UKY- d_Ri*rXR*d_CosTh);
      double ztY =  LC*(kky*UKX - kkx*VKX- d_Ri*rXC*d_CosTh);

      // Construct RHS for density eqn (additional "w" term from background stratification)
      rXR = -LC*(RK[ik].x +(kkz*x0[ik].y-kkx*z0[ik].y)*wki*d_CosTh +(kkx*e0[ik].y-kky*x0[ik].y)*wki*d_SinTh);
      rXC = -LC*(RK[ik].y +(kkx*z0[ik].x-kkz*x0[ik].x)*wki*d_CosTh +(kky*x0[ik].x-kkx*e0[ik].x)*wki*d_SinTh);

      //check if this is the forcing wave number and make appropriate adjustment
      if(kkx==0.0 && kkz==0.0 && kky==d_KFY){ 
	ztX += d_AMPFOR;
      }
 
      // temporary coefficient for Crank-Nicholson step on dissipation terms
      double damp1 = (1.-0.25*d_DELT*d_v2*wk     )/(1.+0.25*d_DELT*d_v2*wk);
      double damp2 =  (1.-0.25*d_DELT*d_v2*wk/d_Sc)/(1.+0.25*d_DELT*d_v2*wk/d_Sc);

      // RK mid step
      x1[ik].x = d_DELT*xiX;
      x1[ik].y = d_DELT*xiY;
      e1[ik].x = d_DELT*etX;
      e1[ik].y = d_DELT*etY;
      z1[ik].x = d_DELT*ztX;
      z1[ik].y = d_DELT*ztY;
      r1[ik].x = d_DELT*rXR;
      r1[ik].y = d_DELT*rXC;

      x0[ik].x = damp1*(xii[ik].x+0.5*d_DELT*xiX);
      x0[ik].y = damp1*(xii[ik].y+0.5*d_DELT*xiY);
      e0[ik].x = damp1*(eta[ik].x+0.5*d_DELT*etX);
      e0[ik].y = damp1*(eta[ik].y+0.5*d_DELT*etY);
      z0[ik].x = damp1*(zet[ik].x+0.5*d_DELT*ztX);
      z0[ik].y = damp1*(zet[ik].y+0.5*d_DELT*ztY);
      r0[ik].x = damp2*(rho[ik].x+0.5*d_DELT*rXR);
      r0[ik].y = damp2*(rho[ik].y+0.5*d_DELT*rXC);

      // Update spectral velocity coeff
      UK[ik].x = -(kky*z0[ik].y-kkz*e0[ik].y)*wki;
      UK[ik].y = -(kkz*e0[ik].x-kky*z0[ik].x)*wki;

      VK[ik].x = -(kkz*x0[ik].y-kkx*z0[ik].y)*wki;
      VK[ik].y = -(kkx*z0[ik].x-kkz*x0[ik].x)*wki;

      WK[ik].x = -(kkx*e0[ik].y-kky*x0[ik].y)*wki;
      WK[ik].y = -(kky*x0[ik].x-kkx*e0[ik].x)*wki;

      ik+=blockDim.x*gridDim.x;
    
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void Step2(cufftDoubleComplex *xii,cufftDoubleComplex *eta,cufftDoubleComplex *zet,cufftDoubleComplex *rho,
		      cufftDoubleComplex *x0,cufftDoubleComplex *e0,cufftDoubleComplex *z0,cufftDoubleComplex *r0,
		      cufftDoubleComplex *x1,cufftDoubleComplex *e1,cufftDoubleComplex *z1,cufftDoubleComplex *r1,
		      cufftDoubleComplex *UK,cufftDoubleComplex *VK,cufftDoubleComplex *WK,cufftDoubleComplex *RK,
		      double *kx,double *ky,double *kz, int *LL){
  int ik = threadIdx.x + blockIdx.x*blockDim.x;
  // This kernel does the first type of sub-step for RK4
  // Most quantities are loaded into local variables to help caching 
  while(ik < d_OK){
      double kky = ky[ik];
      double kkx = kx[ik];
      double kkz = kz[ik];
      double UKX = UK[ik].x;
      double UKY = UK[ik].y;
      double VKX = VK[ik].x;
      double VKY = VK[ik].y;
      double WKX = WK[ik].x;
      double WKY = WK[ik].y;
      double ROX = r0[ik].x;
      double ROY = r0[ik].y;

      double wk =kkx*kkx + kky*kky +kkz*kkz;
      double wki = 1.0/(max(wk,0.001));

      // gradients of density
      double rXR = -kkx*ROY;
      double rXC =  kkx*ROX;
      double rYR = -kky*ROY;
      double rYC =  kky*ROX;
      double rZR = -kkz*ROY;
      double rZC =  kkz*ROX;

      // dealiasing mask (2/3rds or Hou filter, this choice is important for GMRES too!)
      double LC = LL[ik];
      //double LC = exp(-36.0*(pow(2*sqrt(wk)/double(d_NX),36.0)));

      // Do final part of 'convol' (CURL) and add buoyancy term projected due to angle of g
      double xiX =  LC*(kky*WKY - kkz*VKY- d_Ri*(rYR*d_SinTh - rZR*d_CosTh));
      double xiY =  LC*(kkz*VKX - kky*WKX- d_Ri*(rYC*d_SinTh - rZC*d_CosTh));

      double etX =  LC*(kkz*UKY - kkx*WKY+ d_Ri*rXR*d_SinTh);
      double etY =  LC*(kkx*WKX - kkz*UKX+ d_Ri*rXC*d_SinTh);

      double ztX =  LC*(kkx*VKY - kky*UKY- d_Ri*rXR*d_CosTh);
      double ztY =  LC*(kky*UKX - kkx*VKX- d_Ri*rXC*d_CosTh);

      // Construct RHS for density eqn (additional "w" term from background stratification)
      rXR = -LC*(RK[ik].x +(kkz*x0[ik].y-kkx*z0[ik].y)*wki*d_CosTh +(kkx*e0[ik].y-kky*x0[ik].y)*wki*d_SinTh);
      rXC = -LC*(RK[ik].y +(kkx*z0[ik].x-kkz*x0[ik].x)*wki*d_CosTh +(kky*x0[ik].x-kkx*e0[ik].x)*wki*d_SinTh);

      // temporary coefficient for Crank-Nicholson step on dissipation terms
      double damp1 = (1.-0.5*d_DELT*d_v2*wk     )/(1.+0.5*d_DELT*d_v2*wk);
      double damp2 =  (1.-0.5*d_DELT*d_v2*wk/d_Sc)/(1.+0.5*d_DELT*d_v2*wk/d_Sc);

      //check if this is the forcing wave number and make appropriate adjustment
      if(kkx==0.0 && kkz==0.0 && kky==d_KFY){
	ztX += d_AMPFOR;
      }

      // RK mid step
      x1[ik].x = d_DELT*xiX;
      x1[ik].y = d_DELT*xiY;
      e1[ik].x = d_DELT*etX;
      e1[ik].y = d_DELT*etY;
      z1[ik].x = d_DELT*ztX;
      z1[ik].y = d_DELT*ztY;
      r1[ik].x = d_DELT*rXR;
      r1[ik].y = d_DELT*rXC;

      x0[ik].x = damp1*(xii[ik].x+d_DELT*xiX);
      x0[ik].y = damp1*(xii[ik].y+d_DELT*xiY);
      e0[ik].x = damp1*(eta[ik].x+d_DELT*etX);
      e0[ik].y = damp1*(eta[ik].y+d_DELT*etY);
      z0[ik].x = damp1*(zet[ik].x+d_DELT*ztX);
      z0[ik].y = damp1*(zet[ik].y+d_DELT*ztY);
      r0[ik].x = damp2*(rho[ik].x+d_DELT*rXR);
      r0[ik].y = damp2*(rho[ik].y+d_DELT*rXC);

      // Update spectral velocity coeff
      UK[ik].x = -(kky*z0[ik].y-kkz*e0[ik].y)*wki;
      UK[ik].y = -(kkz*e0[ik].x-kky*z0[ik].x)*wki;

      VK[ik].x = -(kkz*x0[ik].y-kkx*z0[ik].y)*wki;
      VK[ik].y = -(kkx*z0[ik].x-kkz*x0[ik].x)*wki;

      WK[ik].x = -(kkx*e0[ik].y-kky*x0[ik].y)*wki;
      WK[ik].y = -(kky*x0[ik].x-kkx*e0[ik].x)*wki;

      ik+=blockDim.x*gridDim.x;
    
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void Step(cufftDoubleComplex *xii,cufftDoubleComplex *eta,cufftDoubleComplex *zet,cufftDoubleComplex *rho,
		     cufftDoubleComplex *x0,cufftDoubleComplex *e0,cufftDoubleComplex *z0,cufftDoubleComplex *r0,
		     cufftDoubleComplex *x1,cufftDoubleComplex *e1,cufftDoubleComplex *z1,cufftDoubleComplex *r1,
		     cufftDoubleComplex *x2,cufftDoubleComplex *e2,cufftDoubleComplex *z2,cufftDoubleComplex *r2,
		     cufftDoubleComplex *x3,cufftDoubleComplex *e3,cufftDoubleComplex *z3,cufftDoubleComplex *r3,
		     cufftDoubleComplex *UK,cufftDoubleComplex *VK,cufftDoubleComplex *WK,cufftDoubleComplex *RK,
		     double *kx,double *ky,double *kz, int *LL){
  int ik = threadIdx.x + blockIdx.x*blockDim.x;
  // Final part of RK4 timestep
  while(ik < d_OK){
      double kky = ky[ik];
      double kkx = kx[ik];
      double kkz = kz[ik];
      double UKX = UK[ik].x;
      double UKY = UK[ik].y;
      double VKX = VK[ik].x;
      double VKY = VK[ik].y;
      double WKX = WK[ik].x;
      double WKY = WK[ik].y;
      double ROX = r0[ik].x;
      double ROY = r0[ik].y;

      double wk  = kkx*kkx + kky*kky +kkz*kkz;
      double wki = 1.0/(max(wk,0.001));

      double rXR = -kkx*ROY;
      double rXC =  kkx*ROX;
      double rYR = -kky*ROY;
      double rYC =  kky*ROX;
      double rZR = -kkz*ROY;
      double rZC =  kkz*ROX;
      
      // temporary coefficient for Crank-Nicholson step on dissipation terms
      double damp1 =  (1.-0.5*d_DELT*d_v2*wk     );
      double damp2 =  (1.-0.5*d_DELT*d_v2*wk/d_Sc);

      // dealiasing mask (2/3rds or Hou filter, this choice is important for GMRES too!)
      double LC = LL[ik];
      //double LC = exp(-36.0*(pow(2*sqrt(wk)/double(d_NX),36.0)));

 
      // Do final part of 'convol' (CURL) and add buoyancy term projected due to angle of g
      double xiX =  kky*WKY - kkz*VKY- d_Ri*(rYR*d_SinTh - rZR*d_CosTh);
      double xiY =  kkz*VKX - kky*WKX- d_Ri*(rYC*d_SinTh - rZC*d_CosTh);

      double etX =  kkz*UKY - kkx*WKY+ d_Ri*rXR*d_SinTh;
      double etY =  kkx*WKX - kkz*UKX+ d_Ri*rXC*d_SinTh;

      double ztX =  kkx*VKY - kky*UKY- d_Ri*rXR*d_CosTh;
      double ztY =  kky*UKX - kkx*VKX- d_Ri*rXC*d_CosTh;

      // Construct RHS for density eqn (additional "w" term from background stratification)
      rXR = -(RK[ik].x +(kkz*x0[ik].y-kkx*z0[ik].y)*wki*d_CosTh +(kkx*e0[ik].y-kky*x0[ik].y)*wki*d_SinTh);
      rXC = -(RK[ik].y +(kkx*z0[ik].x-kkz*x0[ik].x)*wki*d_CosTh +(kky*x0[ik].x-kkx*e0[ik].x)*wki*d_SinTh);
      //check if this is the forcing wave number and make appropriate adjustment

      if(kkx==0.0 && kkz==0.0 && kky==d_KFY){
	ztX += d_AMPFOR;
      }

      // RK final step
      xii[ik].x = LC*(damp1*xii[ik].x + sixth*(x1[ik].x+d_DELT*xiX)+third*(x2[ik].x+x3[ik].x))/(1.+0.5*d_DELT*d_v2*wk);;
      xii[ik].y = LC*(damp1*xii[ik].y + sixth*(x1[ik].y+d_DELT*xiY)+third*(x2[ik].y+x3[ik].y))/(1.+0.5*d_DELT*d_v2*wk);;
      eta[ik].x = LC*(damp1*eta[ik].x + sixth*(e1[ik].x+d_DELT*etX)+third*(e2[ik].x+e3[ik].x))/(1.+0.5*d_DELT*d_v2*wk);;
      eta[ik].y = LC*(damp1*eta[ik].y + sixth*(e1[ik].y+d_DELT*etY)+third*(e2[ik].y+e3[ik].y))/(1.+0.5*d_DELT*d_v2*wk);;
      zet[ik].x = LC*(damp1*zet[ik].x + sixth*(z1[ik].x+d_DELT*ztX)+third*(z2[ik].x+z3[ik].x))/(1.+0.5*d_DELT*d_v2*wk);;
      zet[ik].y = LC*(damp1*zet[ik].y + sixth*(z1[ik].y+d_DELT*ztY)+third*(z2[ik].y+z3[ik].y))/(1.+0.5*d_DELT*d_v2*wk);;
      rho[ik].x = LC*(damp2*rho[ik].x + sixth*(r1[ik].x+d_DELT*rXR)+third*(r2[ik].x+r3[ik].x))/(1.+0.5*d_DELT*d_v2*wk/d_Sc);;
      rho[ik].y = LC*(damp2*rho[ik].y + sixth*(r1[ik].y+d_DELT*rXC)+third*(r2[ik].y+r3[ik].y))/(1.+0.5*d_DELT*d_v2*wk/d_Sc);;

      x0[ik].x = xii[ik].x;
      x0[ik].y = xii[ik].y;
      e0[ik].x = eta[ik].x;
      e0[ik].y = eta[ik].y;
      z0[ik].x = zet[ik].x;
      z0[ik].y = zet[ik].y;
      r0[ik].x = rho[ik].x;
      r0[ik].y = rho[ik].y;

      /* // remove any divergence */
      double phiR =  (kkx*xii[ik].y + kky*eta[ik].y + kkz*zet[ik].y)*wki;
      double phiC = -(kkx*xii[ik].x + kky*eta[ik].x + kkz*zet[ik].x)*wki;

      xii[ik].x -= -kkx*phiC;
      xii[ik].y -=  kkx*phiR;
      eta[ik].x -= -kky*phiC;
      eta[ik].y -=  kky*phiR;
      zet[ik].x -= -kkz*phiC;
      zet[ik].y -=  kkz*phiR;

      // Update spectral velocity coeff
      UK[ik].x = -(kky*z0[ik].y-kkz*e0[ik].y)*wki;
      UK[ik].y = -(kkz*e0[ik].x-kky*z0[ik].x)*wki;

      VK[ik].x = -(kkz*x0[ik].y-kkx*z0[ik].y)*wki;
      VK[ik].y = -(kkx*z0[ik].x-kkz*x0[ik].x)*wki;

      WK[ik].x = -(kkx*e0[ik].y-kky*x0[ik].y)*wki;
      WK[ik].y = -(kky*x0[ik].x-kkx*e0[ik].x)*wki;

      ik+=blockDim.x*gridDim.x;
    
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void RHS(cufftDoubleComplex *xii,cufftDoubleComplex *eta,cufftDoubleComplex *zet,cufftDoubleComplex *rho,
		      cufftDoubleComplex *x0,cufftDoubleComplex *e0,cufftDoubleComplex *z0,cufftDoubleComplex *r0,
		      cufftDoubleComplex *UK,cufftDoubleComplex *VK,cufftDoubleComplex *WK,cufftDoubleComplex *RK,
		      double *kx,double *ky,double *kz, int *LL){
  int ik = threadIdx.x + blockIdx.x*blockDim.x;
  // RHS of prognostic eqns for use in GMRES (i.e. du/dt )
  while(ik < d_OK){
      double kky = ky[ik];
      double kkx = kx[ik];
      double kkz = kz[ik];
      double UKX = UK[ik].x;
      double UKY = UK[ik].y;
      double VKX = VK[ik].x;
      double VKY = VK[ik].y;
      double WKX = WK[ik].x;
      double WKY = WK[ik].y;
      double ROX = r0[ik].x;
      double ROY = r0[ik].y;

      double wk  = kkx*kkx + kky*kky +kkz*kkz;
      double wki = 1.0/(max(wk,0.001));
      
      //compute gradients of density
      double rXR = -kkx*ROY;
      double rXC =  kkx*ROX;
      double rYR = -kky*ROY;
      double rYC =  kky*ROX;
      double rZR = -kkz*ROY;
      double rZC =  kkz*ROX;

      // dealiasing mask (2/3rds or Hou filter, this choice is important for GMRES too!)
      double LC = LL[ik];
      //double LC = exp(-36.0*(pow(2*sqrt(wk)/double(d_NX),36.0)));

      // Do final part of 'convol' (CURL) and add buoyancy term projected due to angle of g
      double xiX =  LC*(kky*WKY - kkz*VKY- d_Ri*(rYR*d_SinTh - rZR*d_CosTh));
      double xiY =  LC*(kkz*VKX - kky*WKX- d_Ri*(rYC*d_SinTh - rZC*d_CosTh));

      double etX =  LC*(kkz*UKY - kkx*WKY+ d_Ri*rXR*d_SinTh);
      double etY =  LC*(kkx*WKX - kkz*UKX+ d_Ri*rXC*d_SinTh);

      double ztX =  LC*(kkx*VKY - kky*UKY- d_Ri*rXR*d_CosTh);
      double ztY =  LC*(kky*UKX - kkx*VKX- d_Ri*rXC*d_CosTh);

      // Construct RHS for density eqn (additional "w" term from background stratification)
      rXR = -LC*(RK[ik].x +(kkz*x0[ik].y-kkx*z0[ik].y)*wki*d_CosTh +(kkx*e0[ik].y-kky*x0[ik].y)*wki*d_SinTh);
      rXC = -LC*(RK[ik].y +(kkx*z0[ik].x-kkz*x0[ik].x)*wki*d_CosTh +(kky*x0[ik].x-kkx*e0[ik].x)*wki*d_SinTh);
      //check if this is the forcing wave number and make appropriate adjustment

      if(kkx==0.0 && kkz==0.0 && kky==d_KFY){ 
	ztX += d_AMPFOR;
      }

      double damp1 = -d_v2*wk;
      double damp2 = -d_v2*wk/d_Sc;

      x0[ik].x = damp1*xii[ik].x+xiX;
      x0[ik].y = damp1*xii[ik].y+xiY;
      e0[ik].x = damp1*eta[ik].x+etX;
      e0[ik].y = damp1*eta[ik].y+etY;
      z0[ik].x = damp1*zet[ik].x+ztX;
      z0[ik].y = damp1*zet[ik].y+ztY;
      r0[ik].x = damp2*rho[ik].x+rXR;
      r0[ik].y = damp2*rho[ik].y+rXC;

      ik+=blockDim.x*gridDim.x;
    
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
