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
	$(MEXDIR)/test_dchud.$(MEXEXT) \
	$(MEXDIR)/test_dchdd.$(MEXEXT)

# test mex executable 
$(MEXDIR)/test_dchud.$(MEXEXT): $(MEXDIR)/test_dchud.o $(SRCDIR)/cholesky.o
	$(CC) $(LDFLAGS) $(LIBS) $(CFLAGS) -o $@ $^
	
$(MEXDIR)/test_dchdd.$(MEXEXT): $(MEXDIR)/test_dchdd.o $(SRCDIR)/cholesky.o
	$(CC) $(LDFLAGS) $(LIBS) $(CFLAGS) -o $@ $^

# mex object files
$(MEXDIR)/test_dchud.o: $(MEXDIR)/test_dchud.cpp $(INCLUDEDIR)/blas_ext.h
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<
	
$(MEXDIR)/test_dchdd.o: $(MEXDIR)/test_dchdd.cpp $(INCLUDEDIR)/blas_ext.h
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
