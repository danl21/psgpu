	GPU Pseudospectral solver for the triply periodic 3D incompressible Euler equations 

	 For referencing please see the following paper:

	"Sustaining processes via recurrent flows in body-forced turbulence"
	 Lucas & Kerswell, J. Fluid Mech. (2017) https://doi.org/10.1017/jfm.2017.97

____________________________________________________________________________________________________________________________

Setting up a job:

	   Resolution and geometry is set at compile time in src/global_params.F90

	   Input files are params.txt for DNS and additional paramsEXTRA.txt for GMRES
	   plus an optional initial condition file, e.g. Zk.in

____________________________________________________________________________________________________________________________

Runtime parameters:

	params.txt

	    ISEED:	random number seed for IC, ensure well separated for each case
	    IREST: 	restart index from file for files with multiple state vectors (if 0 a random IC is generated)
	    infile:	file name/path for input IC
	    infile type: type of input file (A-H) see init.F90 for differences (e.g. UPO or DNS output)
	    TSTEP:	Total simulation duration in non-dimensional units
	    outfreq1:	frequency of high frequency outputs
	    outfreq2:	frequency of low frequency outputs
	    DELT:	timestep
	    AMPV:	rescaling of amplitude of random initial condition 
	    STATSFLAG:  flag to determine wether diagnostics are computed

_______________________________________________________________________________________________________	    
  -----------
  Contacts
  -----------

        Authors: Dan Lucas
        Institution: University of St Andrews
        Email:  dl21@st-andrews.ac.uk
	
October 2024
