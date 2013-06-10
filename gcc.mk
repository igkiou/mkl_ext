MKLDIR = /usr/local/intel/mkl
MATLABDIR = /usr/local/matlab

MATLABARCH = glnxa64
MKLARCH = intel64
MEXEXT = $(shell $(MATLABDIR)/bin/mexext)
MAPFILE = mexFunction.map

MKLLIBS = -L$(MKLDIR)/lib/$(MKLARCH) $(MKLDIR)/lib/$(MKLARCH)/libmkl_solver_ilp64.a -Wl,--start-group -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -Wl,--end-group
MATLABLIBS = -L$(MATLABDIR)/bin/$(MATLABARCH) -lmx -lmex -lmat
RPATH = -Wl,-rpath-link,$(MATLABDIR)/bin/$(MATLABARCH)
LIBS = $(RPATH) $(MATLABLIBS) $(MKLLIBS) -lm -lpthread 
ifdef USE_CUDA
	CULALIBS = -L$(CULADIR)/$(CULAARCH) -lcula_core -lcula_lapack
	LIBS += $(CULALIBS)
	CUDALIBS = -lcublas -lcudart
	LIBS += $(CUDALIBS) 
endif

MKLINCLUDE = -I$(MKLDIR)/include 
MATLABINCLUDE= -I$(MATLABDIR)/extern/include
INCLUDES = $(MKLINCLUDE) $(MATLABINCLUDE)
ifdef USE_CUDA
	CULAINCLUDE= -I$(CULADIR)/include
	INCLUDES += $(CULAINCLUDE)
	CUDAINCLUDE= -I$(CUDADIR)/include
	INCLUDES += $(CUDAINCLUDE)
endif

CC = gcc
MEXFLAGS = -DUSE_MATLAB_INTERFACE -DMATLAB_MEX_FILE -D_GNU_SOURCE -DNDEBUG -fexceptions -fno-omit-frame-pointer
MKLFLAGS = -DMKL_ILP64 -DUSE_DOUBLE_PRECISION -m64
GENERALFLAGS = -fPIC -W -Wall -Wextra -g -pedantic
OPTIMFLAGS = -march=native -O3 -ffast-math -fopenmp -pthread
REPORTSFLAGS = -Winline -Wimplicit
DEBUGFLAG = -g -D__DEBUG__
ifdef DEBUG_MODE
	CFLAGS = $(DEBUGFLAG) $(MEXFLAGS) $(MKLFLAGS) $(GENERALFLAGS) 
else
	CFLAGS = $(MEXFLAGS) $(MKLFLAGS) $(GENERALFLAGS) $(OPTIMFLAGS)
	ifdef PRODUCE_REPORTS
		CFLAGS += $(REPORTSFLAGS) 
	endif
endif 

LDFLAGS = -pthread -shared -Wl,--version-script,$(MATLABDIR)/extern/lib/$(MATLABARCH)/$(MAPFILE) -Wl,--no-undefined
