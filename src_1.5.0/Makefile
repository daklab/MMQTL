# Use system g++ compiler
CC=g++

# Original HPC compiler (commented out)
#CC=/hpc/packages/minerva-common/gcc/7.3.0/bin/g++

DIC=$(PWD)

# Compile flags - using C++11 standard to avoid conflicts
CFLAGS= -c -Wall -g -O3 -std=c++11 -I $(DIC) -I/usr/include/eigen3

# Updated to use system libraries instead of HPC-specific paths
# Using -static for portability across machines
LDFLAGS= -DWITH_LAPACK -static -std=c++11 \
	-I $(DIC)/armadillo-9.400.3/include \
	-I/usr/include/eigen3 \
	-L/usr/lib/gcc/x86_64-linux-gnu/13 \
	-lgsl -lgslcblas \
	-llapack -lblas \
	-lz \
	-lgfortran \
	-fopenmp \
	-lquadmath  

SOURCES=conditional_function.cpp  MeQTLPolyG.cpp  caviar_PostCal.cpp PostCal.cpp  Util.cpp TopKSNP.cpp   gemma.cpp  gemma_gzstream.cpp gemma_io.cpp gemma_lapack.cpp gemma_lmm.cpp gemma_lm.cpp  gemma_mathfunc.cpp gemma_param.cpp  gemma_debug.cpp  gemma_eigenlib.cpp  gemma_fastblas.cpp prepare_input.cpp 

EXECUT=MMQTL26
	
$(EXECUT): $(SOURCES) 
	$(CC) $(SOURCES)   $(LDFLAGS) -o $@
clean:
	rm MMQTL26
