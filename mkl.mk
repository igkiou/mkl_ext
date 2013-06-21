MKLDIR = /usr/local/intel/mkl
MKLARCH = intel64
MKLINT = lp64
ifeq ($(USE_STATIC), 1)
ifeq ($(USE_GCC), 1)
	MKLLIBS = -Wl,--start-group  $(MKLDIR)/lib/$(MKLARCH)/libmkl_intel_$(MKLINT).a $(MKLDIR)/lib/$(MKLARCH)/libmkl_gnu_thread.a $(MKLDIR)/lib/$(MKLARCH)/libmkl_core.a -Wl,--end-group
else
	MKLLIBS = -Wl,--start-group  $(MKLDIR)/lib/$(MKLARCH)/libmkl_intel_$(MKLINT).a $(MKLDIR)/lib/$(MKLARCH)/libmkl_intel_thread.a $(MKLDIR)/lib/$(MKLARCH)/libmkl_core.a -Wl,--end-group
endif
else
ifeq ($(USE_GCC), 1)
	MKLLIBS = -L$(MKLDIR)/lib/$(MKLARCH) -lmkl_intel_$(MKLINT) -lmkl_gnu_thread -lmkl_core
else
	MKLLIBS = -L$(MKLDIR)/lib/$(MKLARCH) -lmkl_intel_$(MKLINT) -lmkl_intel_thread -lmkl_core
endif	
endif
LIBS += $(MKLLIBS)

MKLINCLUDE= -I $(MKLDIR)/include 
INCLUDES += $(MKLINCLUDE)

CFLAGS += -m64
