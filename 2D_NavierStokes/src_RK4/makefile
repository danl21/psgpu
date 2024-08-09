INC     := -I$(CUDA_HOME) -I.
LIB     := -L$(CUDA_HOME)lib64
SDKINC  := -I$(CUDA_HOME)samples/common/inc

NVCCFLAGS := -O3 -arch=sm_70 -D_FORCE_INLINES --use_fast_math -c
#NVCCFLAGS := -O3 -m64 --ptxas-options=-v -arch=sm_20 -c
LIBS    := -lstdc++ -lcufft -lcudart -lm
 
all: cuda dns2 	

cuda:      cuda.cu
	   nvcc cuda.cu $(SDKINC) $(INC) $(NVCCFLAGS) $(LIB) $(LIBS)

dns2: 	global_params.F90 parameters.F90 main.F90 init.F90 cuda.o
	gfortran -O3 -mcmodel=medium -o ../bin/dns2 global_params.F90 parameters.F90 main.F90 init.F90 cuda.o $(LIB) $(LIBS)

rewrite: rewrite_vort rewrite_vort_sym

rewrite_vort:  global_params.F90 rewrite_vort.f 
	gfortran -O3 -mcmodel=large -o ../bin/rewrite_vort global_params.F90 rewrite_vort.f 

rewrite_vort_sym:  global_params.F90 rewrite_vort_sym.f 
	gfortran -O3 -mcmodel=large -o ../bin/rewrite_vort_sym global_params.F90 rewrite_vort_sym.f 

vort_time:  global_params.F90 vort_time.f 
	ifort -O3 -fno-math-errno -mcmodel=large -shared-intel -i8 -r8 -o ../bin/vort_time global_params.F90 vort_time.f 

rewrite_avgs:  global_params.F90 rewrite_avgs.f 
	ifort -O3 -fno-math-errno -mcmodel=large -shared-intel -i8 -r8 -o ../bin/rewrite_avgs global_params.F90 rewrite_avgs.f 

compare_vort:  global_params.F90 compare_vort.f 
	ifort -O3 -fno-math-errno -mcmodel=medium -shared-intel -i8 -r8	 -o ../bin/compare_vort global_params.F90 compare_vort.f 

extract_ZK:  global_params.F90 extract_ZK.f 
	ifort -O3 -fno-math-errno -mcmodel=medium -shared-intel -o ../bin/extract_ZK global_params.F90 extract_ZK.f 

compare_ZK:  global_params.F90 compare_ZK.f 
	ifort -O3 -fno-math-errno -mcmodel=medium -shared-intel -o ../bin/compare_ZK global_params.F90 compare_ZK.f 

rewrite_vort_movie:  global_params.F90 rewrite_vort_movie.f 
	ifort -O3 -fno-math-errno -mcmodel=medium -shared-intel -i8 -r8	 -o ../movie_31/rewrite_vort global_params.F90 rewrite_vort_movie.f

clean:  
	rm -f *.mod *.o

