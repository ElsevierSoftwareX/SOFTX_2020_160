
ifdef _USE_ICC
CPP         = icc
else
CPP         = gcc
endif

UNAME       = $(shell uname -s)
ROOT        = $(ROOTSYS)
ROOTINC     = $(ROOTSYS)/include
WATINC      = ../../wat
LIBODIR     = lib

ifndef HOME_CVODE
   HOME_CVODE       = /Users/valentin/sundials
endif

SRC         = egw_model.cc geodesics.cc merger.cc  numpy.cc GNGen.cc
LIBS        = -lsundials_cvode -lsundials_nvecserial

ifeq ($(UNAME),Darwin)
SHLDFLAGS   = -fno-common -dynamiclib -undefined dynamic_lookup
else
SHLDFLAGS   = -fno-common -shared
endif

CPPSHFLAGS  = 
ifdef _USE_CPP11
CPPSHFLAGS  = -std=c++11
endif
ifdef _USE_ICC
CPPSHFLAGS += -D_USE_ICC
CPPFLAGS   += -D_USE_ICC
CINTOPT    += -D_USE_ICC
endif
ifdef _USE_CPP11
CPPSHFLAGS += -D_USE_CPP11
CPPFLAGS   += -D_USE_CPP11
CINTOPT    += -D_USE_CPP11
endif
ifdef _USE_ROOT6
CPPSHFLAGS += -D_USE_ROOT6
CPPFLAGS   += -D_USE_ROOT6
CINTOPT    += -D_USE_ROOT6
endif
ifdef _USE_LAL
CPPSHFLAGS += -D_USE_LAL
CPPFLAGS   += -D_USE_LAL
CINTOPT    += -D_USE_LAL
endif
ifdef _USE_EBBH
CPPSHFLAGS += -D_USE_EBBH
CPPFLAGS   += -D_USE_EBBH
CINTOPT    += -D_USE_EBBH
endif
ifdef _USE_HEALPIX
CPPSHFLAGS += -D_USE_HEALPIX
CPPFLAGS   += -D_USE_HEALPIX
CINTOPT    += -D_USE_HEALPIX
endif


all:	dir
#	$(CPP)  -I$(HOME_CVODE)/include -L$(HOME_CVODE)/lib  $(LIBS) $(SRC) gw_signal.cc -o gw_signal
	rootcint -f eBBH_dict.cc -c -p -D_USE_ROOT -I$(WATINC) -I$(ROOTINC) GNGen.hh eBBH.hh eBBH_LinkDef.h
	$(CPP) $(SHLDFLAGS) -O2 -fPIC ${CPPSHFLAGS} -I$(HOME_CVODE)/include -I$(WATINC) -I$(ROOTINC) $(SRC) eBBH.cc eBBH_dict.cc -o lib/eBBH.so

dir:
	if [ ! -d $(LIBODIR) ]; then mkdir -p $(LIBODIR) ; fi;
        
clean:
	rm -f lib/*
	rm -f eBBH_dict.cc eBBH_dict.h
	if [ -d $(LIBODIR) ]; then cd $(LIBODIR) ; \rm -f *.so ; fi;
	if [ -d $(LIBODIR) ]; then rmdir $(LIBODIR) ; fi;
        
