CC = icpc
GENERALFLAGS = -fPIC -g -std=c++11
OPTIMFLAGS = -axSSE3,SSE4.1,SSE4.2,AVX -align -O3 -pipe -ipo -fast -parallel -openmp -pthread
WARNFLAGS = -Wall -Wunused-variable -Wcheck -Wextra-tokens -Wformat -Wformat-security -Wmissing-declarations -Wmissing-prototypes -Wpointer-arith -Wreturn-type -Wsign-compare -Wuninitialized
REPORTSFLAGS = -opt-report1 -par-report1 -vec-report1 -openmp-report1
FPFLAGS = -fp-model fast=2 -no-prec-sqrt
#-Winline -Wimplicit
DEBUGFLAG = -g
CFLAGS += $(DEBUGFLAG) $(GENERALFLAGS) $(OPTIMFLAGS) #$(WARNFLAGS)
ifeq ($(DEBUG_MODE), 0)
	CFLAGS += -DNDEBUG
endif
ifeq ($(PRODUCE_REPORTS), 1)
	CFLAGS += $(REPORTSFLAGS) 
endif
LIBS += -lpthread
INCLUDES += -I/usr/include/x86_64-linux-gnu/c++/4.8/
