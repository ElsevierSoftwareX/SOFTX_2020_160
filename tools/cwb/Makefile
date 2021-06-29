#------------------------------------------------------
# Package:     Wavelet Analysis Tool
# File name:   Makefile
#------------------------------------------------------
# This Makefile requires GNU Make.
#------------------------------------------------------
#eval `/home/klimenko/ligotools/bin/use_ligotools`

ifdef _USE_ICC
CPP	    = icc
LD	    = icc
else
CPP	    = gcc
LD	    = gcc
endif
AR	    = ar r

ifndef XIFO
XIFO        = 4
endif

UNAME       = $(shell uname -s)

# Set directory names
#----------------------------
# ROOT system
ROOT        = $(ROOTSYS)
ROOTINC     = $(ROOTSYS)/include
WATINC      = ../../wat
TOOLBOX     = ../toolbox
HISTORY     = ../history
STFT        = ../stft
GWAT        = ../gwat
FRAME       = ../frame
PLUGINS     = plugins

ifdef _USE_HEALPIX
HEALPIX_INC = $(HOME_HEALPIX)/src/cxx/Healpix_cxx
HEALPIX_SUP_INC = $(HOME_HEALPIX)/src/cxx/cxxsupport
endif

ifdef _USE_LAL
LALSUITE_INC  = -I$(LAL_INC) 
LALSUITE_INC += -D__STDC_CONSTANT_MACROS 
#LALSUITE_INC += -DHAVE_CONFIG_H -D__STDC_CONSTANT_MACROS 
endif

ifdef _USE_EBBH
EBBH_INC  = -I../eBBH 
EBBH_INC += -I$(HOME_CVODE)/include
endif

# Data monitoring tool
#DMT         = $(GDS)
# Virgo Frame library
FR          = $(HOME_FRLIB)
# where to place HTML documentation for access via WWW
HTML_DST    = $(HOME)/public_html/FTL

# Recommended option for Ultra SPARC CPU is 'ultrasparc'
#CPU_OPT = -mcpu=ultrasparc
# Recommended option for Alpha CPU is 'ev56' or 'ev6'
#CPU_OPT= -mcpu=ev6
# Recommended option for Intel is 'pentium' or 'pentiumpro'
#CPU_OPT += -mcpu=pentium

CPP_OPT = -Wno-deprecated -Wno-unknown-pragmas -fPIC -mavx

CPPFLAGS    = -O2 -Wall $(CPU_OPT) $(CPP_OPT)
CPPFLAGS   += -m64
#CPPFLAGS    += -funroll-loops

CPPSHFLAGS  = -O2 -fPIC $(CPU_OPT) $(CPP_OPT)
CPPSHFLAGS +=  -fexceptions
ifdef _USE_CPP11
CPPSHFLAGS += -std=c++11
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

SRC_LIST  = cwb2G.cc cwb1G.cc cwb.cc config.cc ced.cc livetime.cc variability.cc wavenoise.cc mdc.cc
#SRC_LIST  = mdc.cc

MISC_SRC =

# extra header files
EXTRA_H    = 

FTL_MACROS = 

ODIR      = obj
SHODIR	  = shared
LIBODIR   = lib

# ------------------------------------------------------
FTL_SRC  = $(SRC_LIST)
FTL_SRC  += $(MISC_SRC)

FTL_H    = $(FTL_SRC:.cc=.hh)

FTL_OBJ  = $(FTL_SRC:%.cc=$(ODIR)/%.o)

FTL_SHO  = $(FTL_SRC:%.cc=$(SHODIR)/%.sho)
FTL_SHO += $(SHODIR)/cwb_dict.sho

WSH_H    += $(FTL_SRC:.cc=.hh)
WSH_H    += xroot.hh
WSH_H    += cwb_LinkDef.h

CPPSHFLAGS += -D_USE_ROOT -I$(ROOTINC) -I$(WATINC) -I$(TOOLBOX) -I$(HISTORY) -I$(STFT) -I$(GWAT) -I$(FRAME) -I$(PLUGINS) -DXIFO=${XIFO} 
CPPFLAGS += -D_USE_ROOT -I$(ROOTINC) -I$(WATINC) -I$(TOOLBOX) -I$(HISTORY) -I$(STFT) -I$(GWAT) -I$(FRAME) -I$(PLUGINS) -DXIFO=${XIFO}
CINTOPT =  -D_USE_ROOT -I$(ROOTINC) -I$(WATINC) -I$(TOOLBOX) -I$(HISTORY) -I$(STFT) -I$(GWAT) -I$(FRAME) -I$(PLUGINS) -DXIFO=${XIFO}
ifdef _USE_ROOT6
CPPSHFLAGS += -D_USE_ROOT6
CPPFLAGS   += -D_USE_ROOT6
CINTOPT    += -D_USE_ROOT6
endif
ifdef _USE_LAL
CPPSHFLAGS += -D_USE_LAL $(LALSUITE_INC)
CPPFLAGS   += -D_USE_LAL $(LALSUITE_INC)
CINTOPT    += -D_USE_LAL $(LALSUITE_INC)
endif
ifdef _USE_EBBH
CPPSHFLAGS += -D_USE_EBBH $(EBBH_INC)
CPPFLAGS   += -D_USE_EBBH $(EBBH_INC)
CINTOPT    += -D_USE_EBBH $(EBBH_INC)
endif

ifdef _USE_HEALPIX
CPPSHFLAGS += -D_USE_HEALPIX -I$(HEALPIX_INC) -I$(HEALPIX_SUP_INC)
CPPFLAGS   += -D_USE_HEALPIX -I$(HEALPIX_INC) -I$(HEALPIX_SUP_INC)
CINTOPT    += -D_USE_HEALPIX -I$(HEALPIX_INC) -I$(HEALPIX_SUP_INC)
endif

LDFLAGS   =
LIBOPT    = -L. -cwb

SHLDFLAGS = -fno-common -dynamiclib -undefined dynamic_lookup

ifeq ($(shell uname), Linux)
SHLDFLAGS = -shared
endif
ifeq ($(UNAME),Darwin)
WDIR		= ${PWD}
SHLDFLAGS+= -fno-common -dynamiclib -undefined dynamic_lookup
else
WDIR		= .
SHLDFLAGS+= -fno-common -dynamiclib
endif

ifdef DMT
CPPFLAGS   += -D_USE_DMT -I$(DMT)/include
CPPSHFLAGS += -D_USE_DMT -I$(DMT)/include
CINTOPT    += -D_USE_DMT -I$(DMT)/include
LIBOPT     += -L$(DMT)/lib -lzlib
endif

ifdef FR
CPPFLAGS  += -D_USE_FR -I$(FR)/src
CPPSHFLAGS  += -D_USE_FR -I$(FR)/src
CINTOPT  += -D_USE_FR -I$(FR)/src 
LIBOPT     += -L$(FR)/$(UNAME) -lFrame
endif

# ------------------------------------------------------

all  : dir xroot xcwb
xcwb : shared 
	chmod +x scripts/*.sh;
	chmod +x scripts/*.csh;

xinfo:
	cd info ; make

script:
	cd scripts ; make

dir:
	if [ ! -d $(ODIR) ]; then mkdir $(ODIR) ; fi;
	if [ ! -d $(SHODIR) ]; then mkdir $(SHODIR) ; fi;
	if [ ! -d $(LIBODIR) ]; then mkdir -p $(LIBODIR) ; fi;
	#if [ ! -d archive ]; then mkdir archive ; fi;
	#if [ ! -d html ]; then mkdir html ; fi;
	#if [ ! -d macro ]; then mkdir macro ; fi;
	#if [ ! -d lib ]; then mkdir lib ; fi;	
	#if [ ! -d include ]; then mkdir macro ; fi;
clean:
	\rm -f lib/cwb.so 
	\rm -f lib/libcwb.a
	if [ -d $(LIBODIR) ]; then cd $(LIBODIR) ; \rm -f *.so ; fi;
	if [ -d $(LIBODIR) ]; then rmdir $(LIBODIR) ; fi;
	if [ -d $(SHODIR) ]; then cd $(SHODIR) ; \rm -f *.sho ; fi;
	if [ -d $(SHODIR) ]; then rmdir $(SHODIR) ; fi;
	if [ -d $(ODIR) ]; then cd $(ODIR); \rm -f *.o ; fi;
	if [ -d $(ODIR) ]; then rmdir $(ODIR) ; fi;
	\rm -f xroot.hh 
	\rm -f cwb_dict.cc cwb_dict.h 

cleanlib:
	\rm -f lib/libcwb.a
	cd $(ODIR); \rm -f *.o

main :  $(ODIR)/main.o $(FTL_OBJ)
	$(LD) $(LDFLAGS) $(ODIR)/main.o -o main $(LIBOPT)

$(ODIR)/main.o :  main.cc $(FTL_H) lib/libcwb.a
	$(CPP) -c $(CPPFLAGS) main.cc -o $(ODIR)/main.o

$(ODIR)/%.o : %.cc $(FTL_H)
	$(CPP) -c $(CPPFLAGS) $< -o $@

$(SHODIR)/%.sho : %.cc $(FTL_H)
	$(CPP) -c $(CPPSHFLAGS) $< -o $@

t%: $(ODIR)/t%.o lib/libcwb.a 
	$(LD) $(LDFLAGS) $< -o $@ $(LIBOPT)

lib: lib/libcwb.a

lib/libcwb.a: $(FTL_OBJ)
	$(AR) $@ $(FTL_OBJ)

shared: lib/cwb.so

lib/cwb.so: $(FTL_SHO)
	\rm -rf lib/cwb.so
	$(LD) $(SHLDFLAGS) -O2 -g -o $(WDIR)/lib/cwb.so $(FTL_SHO)

cwb_dict.cc: $(WSH_H)
	rootcint -f cwb_dict.cc -c -p $(CINTOPT) $(WSH_H)  

# generate HTML documentation
html: html-doc install-html

html-doc: $(FTL_SRC) $(FTL_H)  
	root -b -n -q -l macro/html_doc.C
install:
	if [ ! -d ${prefix}/macro ]; then mkdir ${prefix}/macro ;else \rm -f ${prefix}/macro/* ; fi;
	if [ ! -d ${prefix}/lib ]; then mkdir ${prefix}/lib ;else \rm -f ${prefix}/lib/* ; fi;	
	if [ ! -d ${prefix}/include ]; then mkdir ${prefix}/include ; else \rm -f ${prefix}/include/* ;  fi;
	cp -r macro/*  ${prefix}/macro/ 
	cp -r include/* ${prefix}/include/
	cp -r lib/*  ${prefix}/lib/

# copy html documentation to place where it will be accessed via WWW
install-html:
	if [ -n $(HTML_DST) -a -d $(HTML_DST) ] ; \
        then cp -r html/* $(HTML_DST); fi

# create archives 
tar: archive/cwb.tar archive/macro.tar

archive/wat.tar: $(FTL_SRC) $(MISC_SRC) $(FTL_H) Makefile
	if [ -f $@ ]; then mv $@ $@_old ; fi;
	tar cf $@ $(SRC_LIST) $(MISC_SRC) $(WSH_H) $(EXTRA_H) Makefile INSTALL

archive/macro.tar: macro/*.C
	if [ -f $@ ]; then mv $@ $@_old ; fi;
	tar cf $@ macro/*.C
	chmod a+r $@

# create PostScript Document for compact printout (2 pages on each sheet)

PS_HEADER='$$n|<< Wavelet Analysis Tool -- $$W>>|Page $$% of $$='

ps: 
	enscript --header=$(PS_HEADER) --line-numbers --pretty-print=cpp -o - \
	$(foreach f, $(basename $(SRC_LIST) $(MISC_SRC)),$f.hh $f.cc)\
	$(EXTRA_H) \
	cwb_LinkDef.h Makefile INSTALL | \
        psnup -2 > archive/cwb.ps

xroot:
ifdef _USE_ROOT6
	echo    '#ifndef XROOT_HH'                   >  xroot.hh
	echo    '#define XROOT_HH'                   >> xroot.hh
	echo    ''                                   >> xroot.hh
	echo    '#define SEARCH(TYPE) TYPE cfg_search'  >> xroot.hh
	echo    '#define GAMMA(TYPE)  TYPE cfg_gamma'   >> xroot.hh
	echo    '#define XROOT(ROOT5,ROOT6) ROOT6'   >> xroot.hh
	echo    ''                                   >> xroot.hh
	echo    '#endif'                             >> xroot.hh
else
	echo    '#ifndef XROOT_HH'                   >  xroot.hh
	echo    '#define XROOT_HH'                   >> xroot.hh
	echo    ''                                   >> xroot.hh
	echo    '#define SEARCH(TYPE) TYPE search'   >> xroot.hh
	echo    '#define GAMMA(TYPE)  TYPE gamma'    >> xroot.hh
	echo    '#define XROOT(ROOT5,ROOT6) ROOT5'   >> xroot.hh
	echo    ''                                   >> xroot.hh
	echo    '#endif'                             >> xroot.hh
endif

