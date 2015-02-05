#USE_GCC = 1
CFLAGS =
LDFLAGS =
INCLUDES = -I/home/igkiou/MATLAB/mkl_ext/linalg_benchmarking/armadillo-4.600.4/include -I/home/igkiou/MATLAB/mkl_ext/linalg_benchmarking/blaze-2.2 -I/usr/local/include/boost
LIBS = -L/home/igkiou/MATLAB/mkl_ext/linalg_benchmarking/armadillo-4.600.4/build -larmadillo -lboost_system -lboost_thread

ifeq ($(USE_GCC), 1)
	include gcc.mk
else
	include icc.mk
endif

include matlab.mk
include mkl.mk

CFLAGS += -DARMA_NO_DEBUG #-fwhole-program

all: gemm symm

gemm: gemm_mkl_test.$(MEXEXT) gemm_arma_test.$(MEXEXT) gemm_blaze_test.$(MEXEXT)

symm: symm_mkl_test.$(MEXEXT) symm_arma_test.$(MEXEXT)  symm_blaze_test.$(MEXEXT)

%.$(MEXEXT): %.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) $(INCLUDES) -o $@  $<  $(LIBS) 
		
clean:	
	rm -rf *.o *~

distclean:	
	rm -rf *.o *~ *.$(MEXEXT)
