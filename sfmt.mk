SFMTDIR = SFMT-src-1.4
#ifeq ($(USE_STATIC), 1)
	SFMTLIBS = $(SFMTDIR)/libsfmt.a
#else
#	DSFMTLIBS = -L$(SFMTDIR) -lsfmt
#endif

LIBS += $(SFMTLIBS)

SFMTINCLUDE = -I $(SFMTDIR) 
INCLUDES += $(SFMTINCLUDE)
CFLAGS += -DHAVE_SSE2=1