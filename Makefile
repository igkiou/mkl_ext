# PRODUCE_REPORTS = 1 
# USE_GUIDE = 1
# GENERATE_PROFILE = 1 
# USE_PROFILE = 1
DEBUG_MODE = 1
USE_MKL = 1
# USE_GCC = 1
ifdef USE_GCC
	include gcc.mk
else
	include icc.mk
endif

ROOTDIR = $(shell pwd)
INCLUDEDIR = $(ROOTDIR)/include
LIBDIR = $(ROOTDIR)/lib
SRCDIR = $(ROOTDIR)/src
MEXDIR = $(ROOTDIR)/mexfiles
BINDIR = $(ROOTDIR)/bin
INCLUDES += -I$(INCLUDEDIR)
LIBS += -L$(LIBDIR)

all: tests

#$(LIBDIR)/libmkl_ext.so: $(SRCDIR)/mkl_cholesky.o
#	$(CC) -shared -Wl,-soname,$@.1 -o $@.1.0.0 $^
#	ln -sf  $@.1.0.0 $@.1
#	ln -sf  $@.1 $@

tests: \
	$(MEXDIR)/dchud.$(MEXEXT) \
	$(MEXDIR)/dchdd.$(MEXEXT) \
	$(MEXDIR)/dchr.$(MEXEXT) \
	$(MEXDIR)/dchmv.$(MEXEXT) \
	$(MEXDIR)/dchrk.$(MEXEXT) \
	$(MEXDIR)/dchmm.$(MEXEXT)

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

# mex object files
$(MEXDIR)/dchud.o: $(MEXDIR)/dchud.cpp $(INCLUDEDIR)/blas_ext.h
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<
	
$(MEXDIR)/dchdd.o: $(MEXDIR)/dchdd.cpp $(INCLUDEDIR)/blas_ext.h
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<
	
$(MEXDIR)/dchr.o: $(MEXDIR)/dchr.cpp $(INCLUDEDIR)/blas_ext.h
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<
	
$(MEXDIR)/dchmv.o: $(MEXDIR)/dchmv.cpp $(INCLUDEDIR)/blas_ext.h
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<
	
$(MEXDIR)/dchrk.o: $(MEXDIR)/dchrk.cpp $(INCLUDEDIR)/blas_ext.h
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<
	
$(MEXDIR)/dchmm.o: $(MEXDIR)/dchmm.cpp $(INCLUDEDIR)/blas_ext.h
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
