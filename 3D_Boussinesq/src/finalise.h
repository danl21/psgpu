////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// THIS IS A SOURCE FILE FOR THE 3D_PSPGU CODE                          
// IT CONTAINS SOME WRAPPING UP ON THE HOST SIDE
// INCLUDING COMPUTING AND OUTPUTTING SOME MEANS
// PDFS AND SPECTRA, PLUS FREE-ING ARRAYS
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
  
  Umean=(double*)malloc(sizeof(double)*ny*nz);
  Vmean=(double*)malloc(sizeof(double)*ny*nz);
  Wmean=(double*)malloc(sizeof(double)*ny*nz);
  uhat =(double*)malloc(sizeof(double)*ny*nz);
  vhat =(double*)malloc(sizeof(double)*ny*nz);
  what =(double*)malloc(sizeof(double)*ny*nz);

  if(*statsFLAG==1 && nStop > 1){
    countStats++;
    computeStats(XIIstore,ETAstore,ZETstore,RHOstore,XIItavg,ETAtavg,ZETtavg,RHOtavg,kx,ky,kz,LC,nCopy); //final computation of statistics
    
    // COMPUTE SPECTRA
    // ------------------------------------------------------------------
    double spch[nx];
    double spcv[nx];
    int size=0;
    for(int i=0; i<ikty; i++){
      spch[i]=0.0;
      spcv[i]=0.0;
    }
    for(int i=0; i<nkt; i++){
	double wkh=1e-10;
	double wkv=1e-10;
	int jh=0;
	int jv=0;
	if(*Theta < twopi/8.0){ // do vertical and horizontal spectra depending on orientation
	  wkh = max(sqrt(kx[i]*kx[i] + kz[i]*kz[i]),0.001);
	  wkv = max(ky[i],0.001);
	  jh = int(wkh+0.5);
	  jv = int(wkv+0.5);
	}else{
	  wkh = max(sqrt(kx[i]*kx[i] + ky[i]*ky[i]),0.001);
	  wkv = max(kz[i],0.001);
	  jh = int(wkh+0.5);
	  jv = int(wkv+0.5);
	}
	if(jh > nx || jv > nx){
	  printf("spectrum error\n");
	  continue;
	}
	double wk = max(kx[i]*kx[i] + ky[i]*ky[i] + kz[i]*kz[i],0.001);
	double VZX = XIItavg[i].x*XIItavg[i].x + XIItavg[i].y*XIItavg[i].y;
	double VZY = ETAtavg[i].x*ETAtavg[i].x + ETAtavg[i].y*ETAtavg[i].y;
	double VZZ = ZETtavg[i].x*ZETtavg[i].x + ZETtavg[i].y*ZETtavg[i].y;
	double VZ = (VZX+VZY+VZZ)/wk;
	spch[jh] += VZ;
	spcv[jv] += VZ;
	size = max(size,jh);
	size = max(size,jv);
    }
    for(int j=0; j<size; j++){
      fprintf(spec,"%d %e %e\n",j,spch[j],spcv[j]);
    }
    fflush(spec);
    // ------------------------------------------------------------------

    // WRITE FINAL MEANS ETC.
    // ------------------------------------------------------------------
    // send averages to GPU for processing
    (cudaMemcpy(d_xii,XIItavg,sizeof(cufftDoubleComplex)*nkt,cudaMemcpyHostToDevice));
    (cudaMemcpy(d_eta,ETAtavg,sizeof(cufftDoubleComplex)*nkt,cudaMemcpyHostToDevice));
    (cudaMemcpy(d_zet,ZETtavg,sizeof(cufftDoubleComplex)*nkt,cudaMemcpyHostToDevice));
    (cudaMemcpy(d_rho,RHOtavg,sizeof(cufftDoubleComplex)*nkt,cudaMemcpyHostToDevice));

    // get real space density and write binary file
    KR_FFT(d_rho,d_ro,PlanZ2D,d_ikF,n2);
    (cudaMemcpy(xir,(double *) d_ro, sizeof(double)*nr, cudaMemcpyDeviceToHost));
    fwrite(Time,sizeof(double),1,avgs);
    fwrite(xir,sizeof(double),nr,avgs);

    // get velocity from averaged vorticity and then FFT to real space
    setVelocity<<<nblocks,nthreads>>>(d_xii,d_eta,d_zet,d_rho,d_x0,d_e0,d_z0,d_r0,d_UK,d_VK,d_WK,d_kx,d_ky,d_kz,d_LL);
    KR_FFT_ALL(d_xii,d_eta,d_zet,d_rho,d_UK,d_VK,d_WK,d_xir,d_etr,d_ztr,d_UR,d_VR,d_WR,d_rx,d_ry,d_rz,PlanZ2D,d_ikF,n2);

    // get real space U and write binary file
    (cudaMemcpy(xir,(double *) d_UR, sizeof(double)*nr, cudaMemcpyDeviceToHost));
    fwrite(Time,sizeof(double),1,avgs);
    fwrite(xir,sizeof(double),nr,avgs);

    // get profile of U
    for(int iz=0; iz<nz; iz++){
      for(int iy=0; iy<ny; iy++){
    	int ik = iz*ny+iy;
    	Umean[ik] = 0.0;
	
    	for(int ix=0; ix<nx; ix++){
    	  int ikt = (ny*iz+iy)*nx+ix;
    	  Umean[ik] += xir[ikt];
    	}
    	Umean[ik] /= double(nx);
      }
    }

    // get real space V and write binary file
    (cudaMemcpy(etr,(double *) d_VR, sizeof(double)*nr, cudaMemcpyDeviceToHost));
    fwrite(Time,sizeof(double),1,avgs);
    fwrite(etr,sizeof(double),nr,avgs);

    // get profile of V
    for(int iz=0; iz<nz; iz++){
      for(int iy=0; iy<ny; iy++){
    	int ik = iz*ny+iy;
    	Vmean[ik] = 0.0;
	
    	for(int ix=0; ix<nx; ix++){
    	  int ikt = (ny*iz+iy)*nx+ix;
    	  Vmean[ik] += etr[ikt];
    	}
    	Vmean[ik] /= double(nx);
      }
    }

    // get real space W and write binary file
    (cudaMemcpy(ztr,(double *) d_WR, sizeof(double)*nr, cudaMemcpyDeviceToHost));
    fwrite(Time,sizeof(double),1,avgs);
    fwrite(ztr,sizeof(double),nr,avgs);

    // get profile of W
    for(int iz=0; iz<nz; iz++){
      for(int iy=0; iy<ny; iy++){
        int ik = iz*ny+iy;
        Wmean[ik] = 0.0;

        for(int ix=0; ix<nx; ix++){
          int ikt = (ny*iz+iy)*nx+ix;
          Wmean[ik] += ztr[ikt];
        }
        Wmean[ik] /= double(nx);
      }
    }

    // compute urms 
    subtractReal<<<nblocks,nthreads>>>(d_UR,d_Urms);

    // write binary vrms and profile
    (cudaMemcpy(etr,(double *) d_Vrms, sizeof(double)*nr, cudaMemcpyDeviceToHost));
    fwrite(Time,sizeof(double),1,avgs);
    fwrite(etr,sizeof(double),nr,avgs);
    for(int iz=0; iz<nz; iz++){
      for(int iy=0; iy<ny; iy++){
        int ik = iz*ny+iy;
        vhat[ik] = 0.0;

        for(int ix=0; ix<nx; ix++){
          int ikt = (ny*iz+iy)*nx+ix;
          vhat[ik] += etr[ikt];
        }
        vhat[ik] /= double(nx);
      }
    }

    // write binary wrms and profile
    (cudaMemcpy(etr,(double *) d_Wrms, sizeof(double)*nr, cudaMemcpyDeviceToHost));
    fwrite(Time,sizeof(double),1,avgs);
    fwrite(etr,sizeof(double),nr,avgs);
    for(int iz=0; iz<nz; iz++){
      for(int iy=0; iy<ny; iy++){
        int ik = iz*ny+iy;
        what[ik] = 0.0;

        for(int ix=0; ix<nx; ix++){
          int ikt = (ny*iz+iy)*nx+ix;
          what[ik] += etr[ikt];
        }
        what[ik] /= double(nx);
      }
    }

    // write binary urms and profile
    (cudaMemcpy(ztr,(double *) d_UR, sizeof(double)*nr, cudaMemcpyDeviceToHost));
    fwrite(Time,sizeof(double),1,avgs);
    fwrite(ztr,sizeof(double),nr,avgs);

    for(int iz=0; iz<nz; iz++){
      for(int iy=0; iy<ny; iy++){
        int ik = iz*ny+iy;
        uhat[ik] = 0.0;

        for(int ix=0; ix<nx; ix++){
          int ikt = (ny*iz+iy)*nx+ix;
          uhat[ik] += ztr[ikt];
        }
        uhat[ik] /= double(nx);
    	fprintf(means,"%.15e %.15e %.15e %.15e %.15e %.15e \n",iy*twopi/double(ny),iz*twopi/double(nz),Umean[ik],Vmean[ik],Wmean[ik],uhat[ik],vhat[ik],what[ik]);
      }
      fprintf(means," \n");
    }

  printf("AVERAGE STATS *****\n");
  printf("energy = %e,  eturb = %e \n",enerAvg,eturbAvg);
  printf("dissipation = %e,  dturb = %e \n",dissAvg,dturbAvg);
  printf("input E = %e \n",inputAvg);
  printf("bouyancy flux = %e \n",bfluxAvg);

  printf("max/min dissipation %e, %e \n ", dissMax,dissMin);
  printf("max/min energy %e, %e \n", enerMax,enerMin);
  printf("max/min input %e, %e \n", inputMax,inputMin);

  // write sweep for DNS processing of continuation run
  fprintf(sweep,"%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e \n",1.0/v2,*Ri,*alpha,dissAvg,inputAvg,mixAvg,RiGMax,pedAvg);
  fflush(sweep);

  // compute and write the pdfs of D,I,E
  // ------------------------------------------------------------------
  if(dissMax-dissMin > 0.001 &&  inputMax-inputMin> 0.001 && enerMax-enerMin>0.001){
    int nBin = 100;
    double correct = 1.0000001;
    int I_E[nBin],I_D[nBin],I_I[nBin];

    double dD = (dissMax*correct-dissMin)/double(nBin);
    double dI = (inputMax*correct-inputMin)/double(nBin);
    double dE = (enerMax*correct-enerMin)/double(nBin);

    for(int ib=0; ib<nBin; ib++){
      I_E[ib] = 0;
      I_D[ib] = 0;
      I_I[ib] = 0;
    }
    printf("dD = %e \n",dD);
    printf("dI = %e \n",dI);
    printf("dE = %e \n",dE);
    int j=0;
    for(int NT = 0; NT < countStats; NT++){
      
      j=int((Energy[NT]-enerMin)/dE) +1;
      I_E[j]++;
      j=int((Diss[NT]-dissMin)/dD) +1;
      if(j>nBin) printf("error in j!");
      I_D[j]++;
      j=int((Input[NT]-inputMin)/dI) +1;
      if(j>nBin) printf("error in j!");
      I_I[j]++;
    }
    fflush(stdout);

    for(int ib=0; ib<nBin; ib++){
      fprintf(pdfs,"%.15e %.15e %.15e %.15e %.15e %.15e \n", enerMin+(ib-0.5)*dE , I_E[ib]/((nStats-1)*dE)
  	      , dissMin+(ib-0.5)*dD , I_D[ib]/((nStats-1)*dD)
  	      ,inputMin+(ib-0.5)*dI , I_I[ib]/((nStats-1)*dI));
    }
  }
  }

  // free all of the allocated memory
  // ------------------------------------------------------------------
  free(Umean);
  free(Vmean);
  free(Wmean);
  free(uhat);
  free(vhat);
  free(what);

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
  //Free CPU memory
  if(*statsFLAG !=0){
    free(XIIstore);
    free(ETAstore);
    free(ZETstore);
    free(RHOstore);

    free(XIItavg);
    free(ETAtavg);
    free(ZETtavg);
    free(RHOtavg);

    free(Diss);
    free(Input);
    free(Energy);

    (cudaFree(d_Urms));
    (cudaFree(d_Vrms));
    (cudaFree(d_Wrms));
  
    fclose(sweep);
    fclose(RiG);
    fclose(stats);
    fclose(spec);
    fclose(avgs);
    fclose(means);
    fclose(pdfs);
  }
  printf("stats free done \n"); 
  //Destroy fft plans
  (cufftDestroy(PlanZ2D));
  (cufftDestroy(PlanD2Z));
   
  printf("time stepping done \n"); 
  fflush(stdout);
