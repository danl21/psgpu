____________________________________________________________________________________________________________________________


		GPU Pseudospectral solver for the doubly periodic 2D Navier-Stokes Equations

		For referencing please see the following papers:

                "Spatiotemporal dynamics in 2D Kolmogorov flow over large domains"
		                Lucas & Kerswell, J. Fluid Mech. (2014)

____________________________________________________________________________________________________________________________


Outline of the code:

	   Code is of a standard pseudospectral form where timestepping occurs on fourier coefficients of vorticity
	   and density with nonlinear convolutions computed in physical space. Timestepping scheme is RK4
	   (fourth order Runge-Kutta) with Crank-Nicolson on diffusion terms. 
	   
	   Code is for use on NVIDIA GPUs via the CUDA programming API.
	   
	   Use is made of the CUFFT library for the Fast Fourier Transform (FFT) on GPUs. We also
	   make use of the Thrust library of CUDA kernels for global reductions (spatial averages
	   and maxima). 

	   All DNS routines are found in src/ and include various diagnostic functions and
	   recurrence searching routines for finding close cycles of dynamics for convergence to unstable
	   periodic orbits. 

	   Also included is the Newton-GMRES-hookstep code for converging nonlinear unstable solutions, 
	   inluding periodic orbits. See src/src_gmres for the source code. 

	   Makefiles are included but certain environment variables are required to link to the libs
____________________________________________________________________________________________________________________________

bin/ contains input data, scripts and a short readme for setting up a job.
Please feel free to contact me for assistance or advice on using the code.

  -----------
  Contact
  -----------

        Authors: Dan Lucas
        Institution: Keele University
        Email:  d.lucas1@keele.ac.uk

May 2019
