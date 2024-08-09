INC     := -I$(CUDA_HOME) -I.
LIB     := -L$(CUDA_HOME)lib64
SDKINC  := -I$(CUDA_HOME)samples/common/inc
FFTLIB  := /usr/local/lib

NVCCFLAGS := -O3 -arch=sm_70 -D_FORCE_INLINES --use_fast_math -c
LIBS    := -lstdc++ -lcufft -lcudart -lm

all: cuda gmres

cuda: cuda.cu
	nvcc cuda.cu $(SDKINC) $(INC) $(NVCCFLAGS) $(LIB) $(LIBS)

gmres: 	global_params.F90 parameters.F90 io_arpack.F90 main.F90 timestep.F90 timederiv.F90 convol.F90 fft-fftw.F90 init.F90 out.F90 probe.F90 spec.F90 transf.F90
	gfortran -O3 -mcmodel=medium -o ../../bin/gmres global_params.F90 parameters.F90 io_arpack.F90 main.F90 timestep.F90 timederiv.F90 convol.F90 fft-fftw.F90 init.F90 out.F90 probe.F90 spec.F90 transf.F90 cuda.o -L/home/dan/libs/lapack-3.9.0 -llapack -lrefblas -L/home/dan/libs/ARPACK -larpack $(LIB) $(LIBS) -L$(FFTLIB) -lfftw3 

clean: 
	rm -f *.mod *.o

