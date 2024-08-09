// ------------------------------------------------------------------
// Low frequency outputs: RiG, full 3D field and updates of dt etc.
// ------------------------------------------------------------------

    double tcrit=Tout*double(int((tme+small)/Tout));//want outputs every Tout!                                        
    double tcheck=(tme-delt-tcrit+small)*(tme+small-tcrit);
    
    if(tcheck < 0.0 && *statsFLAG !=0){
      avgVelocity<<<nblocks,nthreads>>>(d_UR,d_VR,d_WR,d_Urms,d_Vrms,d_Wrms,tavgCount); // average real space velocity

      if(*adaptFLAG==1){ // update timestep
	Umax = max(*(thrust::max_element(d_Uptr,d_Uptr+nr)),fabs(*(thrust::min_element(d_Uptr,d_Uptr+nr))));
        delt = max(C/abs(Umax),deltp);
        (cudaMemcpyToSymbol(d_DELT,&delt,sizeof(double)));
      }
    }

    if(tcheck < 0.0 && *statsFLAG !=0 && timestep > 0){
// ------------------------------------------------------------------
      // compute ZT-plane gradient Richarson number (depends on theta)
// ------------------------------------------------------------------
      if(*Theta < twopi/8.0){
	KR_FFT(d_rho,d_ro,PlanZ2D,d_ikF,n2);
	(cudaMemcpy(ztr, d_WR, sizeof(double)*nr, cudaMemcpyDeviceToHost));
	(cudaMemcpy(xir, d_UR, sizeof(double)*nr, cudaMemcpyDeviceToHost));
	(cudaMemcpy(ro, d_ro, sizeof(double)*nr, cudaMemcpyDeviceToHost));
	RiGMax=0.0;
	for(int iy=0; iy<ny; iy++){
	  double Navg=0.0;
	  double Savg=0.0;
	  int iyp = iy+1;
	  if(iy==ny-1) iyp=1;
	  int iym = iy-1;
	  if(iy==0) iyp=ny-1;
	  for(int ix=0; ix<nx; ix++){
	    for(int iz=0; iz<nz; iz++){
	      int ikp = (iz*ny+iyp)*nx+ix;
	      int ikm = (iz*ny+iym)*nx+ix;
	      Navg += 0.5*ny*(ro[ikp]-ro[ikm])/twopi +1.0;
	      Savg += pow(0.5*ny*(xir[ikp]-xir[ikm])/twopi,2);//+ny*(ztr[ikpm]-ztr[ikm])/twopi ;
	    }
	  }
	  fprintf(RiG,"%e %d %e %e\n",tme,iy,Navg/double(nx*nz),Savg/double(nx*nz));
	  RiGMax = max(RiGMax, Navg/double(nx)/Savg/double(nx));
	}
	fprintf(RiG,"   \n");
	fflush(RiG);
      }else{
	KR_FFT(d_rho,d_ro,PlanZ2D,d_ikF,n2);
	(cudaMemcpy(ztr, d_VR, sizeof(double)*nr, cudaMemcpyDeviceToHost));
	(cudaMemcpy(xir, d_UR, sizeof(double)*nr, cudaMemcpyDeviceToHost));
	(cudaMemcpy(ro, d_ro, sizeof(double)*nr, cudaMemcpyDeviceToHost));
	RiGMax=0.0;

	for(int iz=0; iz<nz; iz++){
	  double Navg=0.0;
	  double Savg=0.0;
	  int izp = iz+1;
	  if(iz==nz-1) izp=0;
	  int izm = iz-1;
	  if(iz==0) izm=nz-1;
	  for(int ix=0; ix<nx; ix++){
	    for(int iy=0; iy<ny; iy++){
	      int ikp = (izp*ny+iy)*nx+ix;
	      int ikm = (izm*ny+iy)*nx+ix;
	      Navg += 0.5*nz*(ro[ikp]-ro[ikm])/twopi +1.0;
	      Savg += pow(0.5*nz*(ztr[ikp]-ztr[ikm])/twopi,2) + pow(0.5*nz*(xir[ikp]-xir[ikm])/twopi,2);
	    }
	  }
	  fprintf(RiG,"%e %d %e %e\n",tme,iz,Navg/double(nx*nz),Savg/double(nx*nz));
          RiGMax = max(RiGMax, Navg/double(nx)/Savg/double(nx));
	}
	fprintf(RiG,"   \n");
	fflush(RiG);
      }
    }

// ------------------------------------------------------------------
    double vcrit=Vout*double(int((tme+small)/Vout));//want outputs every Tout!                                        
    if((tme-delt-vcrit+small)*(tme+small-vcrit) < 0.0  && nStop > 1 && *statsFLAG !=0){ 
// ------------------------------------------------------------------
      // output full 3D components
// ------------------------------------------------------------------
      printf("DT = %e \n",delt);
      KR_FFT(d_rho,d_ro,PlanZ2D,d_ikF,n2);
      count++;
      sprintf(xiifile, "JOB%1d/U%03d.raw",rank,count);
      sprintf(etafile, "JOB%1d/V%03d.raw",rank,count);
      sprintf(zetfile, "JOB%1d/W%03d.raw",rank,count);
      sprintf(rhofile, "JOB%1d/rho%03d.raw",rank,count);
      xiiStream = fopen(xiifile,"wb");
      etaStream = fopen(etafile,"wb");
      zetStream = fopen(zetfile,"wb");
      rhoStream = fopen(rhofile,"wb");
      (cudaMemcpy(ztr, d_WR, sizeof(double)*nr, cudaMemcpyDeviceToHost));
      (cudaMemcpy(xir, d_UR, sizeof(double)*nr, cudaMemcpyDeviceToHost));
      (cudaMemcpy(etr, d_VR, sizeof(double)*nr, cudaMemcpyDeviceToHost));
      (cudaMemcpy(ro, d_ro, sizeof(double)*nr, cudaMemcpyDeviceToHost));

      fwrite(xir,sizeof(double),nr,xiiStream);
      fwrite(etr,sizeof(double),nr,etaStream);
      fwrite(ztr,sizeof(double),nr,zetStream);
      fwrite(ro,sizeof(double),nr,rhoStream);

      fclose(etaStream);
      fflush(etaStream);
      fclose(zetStream);
      fflush(zetStream);
      fclose(xiiStream);
      fflush(xiiStream);
      fclose(rhoStream);
      fflush(rhoStream);

      // do components of buoyancy flux
      char bfluxfile[sizeof "JOB#/Bflux000.dat"];
      sprintf(bfluxfile, "JOB%1d/bflux%03d.dat",rank,count);
      bfile = fopen(bfluxfile,"w");

      computeBflux(XIIstore,ETAstore,ZETstore,RHOstore,kx,ky,kz,LC,nCopy);      
    }

// ------------------------------------------------------------------
// ------------------------------------------------------------------
// ------------------------------------------------------------------
