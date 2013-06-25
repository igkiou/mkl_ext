# USE_GCC = 1
USE_BLAS = MKL
USE_RNG = VSL
USE_RNG_NORMAL = BOXMULLER

ROOTDIR = $(shell pwd)
INCLUDEDIR = $(ROOTDIR)/include
LIBDIR = $(ROOTDIR)/lib
SRCDIR = $(ROOTDIR)/src
MEXDIR = $(ROOTDIR)/mexfiles
BINDIR = $(ROOTDIR)/bin
INCLUDES = -I$(INCLUDEDIR)
LIBS = -L$(LIBDIR)

ifeq ($(USE_GCC), 1)
	include gcc.mk
else
	include icc.mk
endif

include matlab.mk

ifeq ($(USE_BLAS), MKL)
	include mkl.mk
	CFLAGS += -DUSE_BLAS_MKL
endif

ifeq ($(USE_RNG), VSL)
	include vsl.mk
	CFLAGS += -DUSE_RNG_VSL
else ifeq ($(USE_RNG), SFMT)
	include sfmt.mk
	CFLAGS += -DUSE_RNG_SFMT
else ifeq ($(USE_RNG), DSFMT)
	include dsfmt.mk
	CFLAGS += -DUSE_RNG_DSFMT
endif

ifeq ($(USE_RNG_NORMAL), MARSAGLIA)
	CFLAGS += -DUSE_RNG_MARSAGLIA
else ifeq ($(USE_RNG_NORMAL), BOXMULLER)
	CFLAGS += -DUSE_RNG_BOX_MULLER
endif

all: tests

#$(LIBDIR)/libmkl_ext.so: $(SRCDIR)/mkl_cholesky.o
#	$(CC) -shared -Wl,-soname,$@.1 -o $@.1.0.0 $^
#	ln -sf  $@.1.0.0 $@.1
#	ln -sf  $@.1 $@

tests: tests_cholesky tests_rng
	
tests_cholesky: \
	$(MEXDIR)/dchud.$(MEXEXT) \
	$(MEXDIR)/dchdd.$(MEXEXT) \
	$(MEXDIR)/dchr.$(MEXEXT) \
	$(MEXDIR)/dchmv.$(MEXEXT) \
	$(MEXDIR)/dchrk.$(MEXEXT) \
	$(MEXDIR)/dchmm.$(MEXEXT) \
	$(MEXDIR)/dchex.$(MEXEXT)

tests_rng: \
	$(MEXDIR)/drnorm.$(MEXEXT) \
	$(MEXDIR)/drunif.$(MEXEXT)

# test mex executable 
$(MEXDIR)/dchud.$(MEXEXT): $(MEXDIR)/dchud.o $(SRCDIR)/cholesky.o
	$(CC) $(LDFLAGS) $(LIBS) $(CFLAGS) -o $@ $^
	
$(MEXDIR)/dchdd.$(MEXEXT): $(MEXDIR)/dchdd.o $(SRCDIR)/cholesky.o
	$(CC) $(LDFLAGS) $(LIBS) $(CFLAGS) -o $@ $^
	
$(MEXDIR)/dchr.$(MEXEXT): $(MEXDIR)/dchr.o $(SRCDIR)/cholesky.o
	$(CC) $(LDFLAGS) $(LIBS) $(CFLAGS) -o $@ $^
	
$(MEXDIR)/dchmv.$(MEXEXT): $(MEXDIR)/dchmv.o $(SRCDIR)/cholesky.o
	$(CC) $(LDFLAGS) $(LIBS) $(CFLAGS) -o $@ $^
	
$(MEXDIR)/dchrk.$(MEXEXT): $(MEXDIR)/dchrk.o $(SRCDIR)/cholesky.o
	$(CC) $(LDFLAGS) $(LIBS) $(CFLAGS) -o $@ $^
	
$(MEXDIR)/dchmm.$(MEXEXT): $(MEXDIR)/dchmm.o $(SRCDIR)/cholesky.o
	$(CC) $(LDFLAGS) $(LIBS) $(CFLAGS) -o $@ $^
	
$(MEXDIR)/dchex.$(MEXEXT): $(MEXDIR)/dchex.o $(SRCDIR)/cholesky.o
	$(CC) $(LDFLAGS) $(LIBS) $(CFLAGS) -o $@ $^
	
$(MEXDIR)/drnorm.$(MEXEXT): $(MEXDIR)/drnorm.o $(SRCDIR)/rng.o
	$(CC) $(LDFLAGS) $(LIBS) $(CFLAGS) -o $@ $^
	
$(MEXDIR)/drunif.$(MEXEXT): $(MEXDIR)/drunif.o $(SRCDIR)/rng.o
	$(CC) $(LDFLAGS) $(LIBS) $(CFLAGS) -o $@ $^
	
# mex object files
$(MEXDIR)/%.o: $(MEXDIR)/%.cpp $(INCLUDEDIR)/blas_ext.h
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<
	
# src object files
$(SRCDIR)/%.o: $(SRCDIR)/%.cpp $(INCLUDEDIR)/%.h $(INCLUDEDIR)/blas_ext.h
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<
	
clean:	
	rm -rf *.o *~
	rm -rf $(MEXDIR)/*.o $(MEXDIR)/*~
	rm -rf $(SRCDIR)/*.o $(SRCDIR)/*~
	rm -rf $(LIBDIR)/*.o $(LIBDIR)/*~

distclean:	
	rm -rf *.o *~
	rm -rf $(MEXDIR)/*.o $(MEXDIR)/*.$(MEXEXT) $(MEXDIR)/*~
	rm -rf $(SRCDIR)/*.o $(SRCDIR)/*~
	rm -rf $(LIBDIR)/*.o $(LIBDIR)/*~ $(LIBDIR)/*.so*
