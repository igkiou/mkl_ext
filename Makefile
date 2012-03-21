# PRODUCE_REPORTS = 1 
# USE_GUIDE = 1
# GENERATE_PROFILE = 1 
# USE_PROFILE = 1
DEBUG_MODE = 1

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

all: mkl_ext tests
	
tests: $(MEXDIR)/test_cholesky.$(MEXEXT)

mkl_ext: $(SRCDIR)/mkl_cholesky.o
	$(CC) -shared -Wl,-soname,$(LIBDIR)/lib$@.so.1 -o $(LIBDIR)/lib$@.so.1.0.0 $^
	ln -sf  $(LIBDIR)/lib$@.so.1.0.0 $(LIBDIR)/lib$@.so.1
	ln -sf  $(LIBDIR)/lib$@.so.1 $(LIBDIR)/lib$@.so
	
## test mex executable 
$(MEXDIR)/test_cholesky.$(MEXEXT): $(MEXDIR)/test_cholesky.o 
	$(CC) $(LDFLAGS) $(LIBS) $(CFLAGS) -lmkl_ext -o $@ $<
#
# mex object files
$(MEXDIR)/%.o: $(MEXDIR)/%.cpp $(INCLUDEDIR)/mkl_cholesky.h
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<
	
# src object files
$(SRCDIR)/%.o: $(SRCDIR)/%.cpp $(INCLUDEDIR)/%.h
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
