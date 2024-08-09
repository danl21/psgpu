INC     := -I$(CUDA_PATH)  -I.
LIB     := -L$(CUDA_PATH)lib64 
SDKINC  := -I$(CUDA_PATH)samples/common/inc

NVCCFLAGS := -O3 -arch=sm_20 --use_fast_math -c
LIBS    := -lcufft -lcudart -lm -lstdc++
  
all: dns rewrite

dns: cuda dns3

cuda:      cuda.cu kernels.h functions.h
	   nvcc cuda.cu $(SDKINC) $(INC) $(NVCCFLAGS) $(LIB) $(LIBS)

dns3: 	global_params.F90 parameters.F90 main.F90 init.F90 cuda.o
	mpif90 -O3 -fno-math-errno -mcmodel=large  -o ../bin/dns_bous global_params.F90 parameters.F90 main.F90 init.F90 cuda.o $(LIB) $(LIBS)

rewrite: rewrite_plane rewrite_avgs

rewrite_plane:  global_params.F90 rewrite_plane.f 
	gfortran -O3 -fno-math-errno -mcmodel=large    -o ../bin/rewrite_plane global_params.F90 rewrite_plane.f 

rewrite_avgs:  global_params.F90 rewrite_avgs.f 
	gfortran -O3 -fno-math-errno -mcmodel=large    -o ../bin/rewrite_avgs global_params.F90 rewrite_avgs.f 

compute_RiG:  global_params.F90 compute_RiG.f 
	gfortran -O3 -fno-math-errno -mcmodel=large    -o ../bin/compute_RiG global_params.F90 compute_RiG.f 

RiG_time:  global_params.F90 RiG_time.f 
	gfortran -O3 -fno-math-errno -mcmodel=large    -o ../bin/RiG_time global_params.F90 RiG_time.f 

RiG_avgs:  global_params.F90 RiG_avgs.f 
	gfortran -O3 -fno-math-errno -mcmodel=large    -o ../bin/RiG_avgs global_params.F90 RiG_avgs.f 

clean: 
	rm *.o ../bin/dns_* ../bin/*write* 
