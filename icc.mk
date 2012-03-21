MKLDIR = /usr/local/intel/mkl
MATLABDIR = /usr/local/matlabR2010a

MATLABARCH = glnxa64
MKLARCH = intel64
MEXEXT = $(shell $(MATLABDIR)/bin/mexext)
MAPFILE = mexFunction.map

MKLLIBS = -L$(MKLDIR)/lib/$(MKLARCH) $(MKLDIR)/lib/$(MKLARCH)/libmkl_solver_ilp64.a -Wl,--start-group -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group
MATLABLIBS = -L$(MATLABDIR)/bin/$(MATLABARCH) -lmx -lmex -lmat
RPATH = -Wl,-rpath-link,$(MATLABDIR)/bin/$(MATLABARCH)
LIBS = $(RPATH) $(MATLABLIBS) $(MKLLIBS) -lm -lpthread

MKLINCLUDE = -I$(MKLDIR)/include 
MATLABINCLUDE= -I$(MATLABDIR)/extern/include
INCLUDES = $(MKLINCLUDE) $(MATLABINCLUDE)

CC = icc
MEXFLAGS = -DUSE_MATLAB_INTERFACE -DMATLAB_MEX_FILE -D_GNU_SOURCE -DNDEBUG -fexceptions -fno-omit-frame-pointer
MKLFLAGS = -DMKL_ILP64 -DUSE_DOUBLE_PRECISION
GENERALFLAGS = -fPIC -g -Wall -Wunused-variable -Wcheck -Wextra-tokens -Wformat -Wformat-security -Wmissing-declarations -Wmissing-prototypes -Wpointer-arith -Wreturn-type -Wsign-compare -Wuninitialized
OPTIMFLAGS = -axSSE3,SSE4.1,SSE4.2,AVX -O3 -pipe -ipo -fast -parallel -openmp -pthread
REPORTSFLAGS = -opt-report 3 -openmp-report2 -par-report3 -vec-report3 -Winline -Wimplicit
FPFLAGS = -fp-model fast=2 -no-prec-sqrt
GUIDEFLAG = -guide=4
PROFGENFLAG = -prof-gen -profile-functions -profile-loops
PROFUSEFLAG = -prof-use
DEBUGFLAG = -g -D__DEBUG__
ifdef DEBUG_MODE
	CFLAGS = $(DEBUGFLAG) $(MEXFLAGS) $(MKLFLAGS) $(GENERALFLAGS) -pthread -openmp
else
	CFLAGS = $(MEXFLAGS) $(MKLFLAGS) $(GENERALFLAGS) $(OPTIMFLAGS)
	ifdef PRODUCE_REPORTS
		CFLAGS += $(REPORTSFLAGS) 
	endif
	ifdef USE_GUIDE
		CFLAGS += $(GUIDEFLAG) 
	endif
	ifdef GENERATE_PROFILE
		CFLAGS += $(PROFGENFLAG) 
	endif
	ifdef USE_PROFILE
		CFLAGS += $(PROFUSEFLAG) 
	endif
endif

LDFLAGS = -pthread -shared -Wl,--version-script,$(MATLABDIR)/extern/lib/$(MATLABARCH)/$(MAPFILE) -Wl,--no-undefined
