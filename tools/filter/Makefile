#------------------------------------------------------
# Package:     Filter Tool
# File name:   Makefile
#------------------------------------------------------
# This Makefile requires GNU Make.
#------------------------------------------------------

ifdef _USE_ICC
CPP         = icc
LD          = icc
else
CPP         = gcc
LD          = gcc
endif
AR	    = ar r

UNAME       = $(shell uname -s)

# Set directory names
#----------------------------
# ROOT system
ROOT        = $(ROOTSYS)
ROOTINC     = $(ROOTSYS)/include
WATINC      = ../../wat

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

CPP_OPT = -Wno-deprecated -fPIC

CPPFLAGS    = -O2 -Wall $(CPU_OPT) $(CPP_OPT)
#CPPFLAGS    += -funroll-loops

CPPSHFLAGS  = -O2 -fPIC $(CPU_OPT) $(CPP_OPT)
CPPSHFLAGS +=  -fexceptions
ifdef _USE_CPP11
CPPSHFLAGS += -std=c++11
CPPFLAGS += -std=c++11
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

SRC_LIST  =	Filter.cc

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
FTL_OBJ += $(ODIR)/mkfcomplex.o

FTL_SHO  = $(FTL_SRC:%.cc=$(SHODIR)/%.sho)
FTL_SHO += $(SHODIR)/Filter_dict.sho
FTL_SHO += $(SHODIR)/mkfcomplex.sho

WSH_H    += $(FTL_SRC:.cc=.hh)
WSH_H    += Filter_LinkDef.h

CPPSHFLAGS += -D_USE_ROOT -I$(ROOTINC) -I$(WATINC)
CPPFLAGS += -D_USE_ROOT -I$(ROOTINC) -I$(WATINC)
CINTOPT =  -D_USE_ROOT -I$(ROOTINC) -I$(WATINC)
ifdef _USE_WAT
CPPSHFLAGS  += -D_USE_WAT
CPPFLAGS    += -D_USE_WAT
CINTOPT     += -D_USE_WAT
endif
LDFLAGS   =
LIBOPT    = -L. -lFilter
SHLDFLAGS = -shared
ifeq ($(UNAME),Darwin)
WDIR	=	${PWD}
SHLDFLAGS+= -fno-common -dynamiclib -undefined dynamic_lookup
else
WDIR	=	.
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
#FTL_SHO += $(FR)/lib/libFrame.so
endif
# ------------------------------------------------------

all:
	make clean; make libs _USE_WAT=1; make rename _USE_WAT=1; make cleanobj; make libs
 
libs: dir lib shared 

rename:
ifdef _USE_WAT
	mv lib/libFilter.a lib/libWFilter.a
	mv lib/Filter.so lib/WFilter.so
endif

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
	\rm -f lib/WFilter.so 
	\rm -f lib/libWFilter.a
	\rm -f lib/Filter.so 
	\rm -f lib/libFilter.a
	if [ -d $(LIBODIR) ]; then cd $(LIBODIR) ; \rm -f *.so ; fi;
	if [ -d $(LIBODIR) ]; then rmdir $(LIBODIR) ; fi;
	if [ -d $(SHODIR) ]; then cd $(SHODIR) ; \rm -f *.sho ; fi;
	if [ -d $(SHODIR) ]; then rmdir $(SHODIR) ; fi;
	if [ -d $(ODIR) ]; then cd $(ODIR); \rm -f *.o ; fi;
	if [ -d $(ODIR) ]; then rmdir $(ODIR) ; fi;
	\rm -f Filter_dict.cc Filter_dict.h 

cleanobj:
	if [ -d $(SHODIR) ]; then cd $(SHODIR) ; \rm -f *.sho ; fi;
	if [ -d $(SHODIR) ]; then rmdir $(SHODIR) ; fi;
	if [ -d $(ODIR) ]; then cd $(ODIR); \rm -f *.o ; fi;
	if [ -d $(ODIR) ]; then rmdir $(ODIR) ; fi;
	\rm -f Filter_dict.cc Filter_dict.h 

cleanlib:
	\rm -f lib/libFilter.a
	cd $(ODIR); \rm -f *.o

main :  $(ODIR)/main.o $(FTL_OBJ)
	$(LD) $(LDFLAGS) $(ODIR)/main.o -o main $(LIBOPT)

$(ODIR)/main.o :  main.cc $(FTL_H) lib/libFilter.a
	$(CPP) -c $(CPPFLAGS) main.cc -o $(ODIR)/main.o

$(ODIR)/%.o : %.cc $(FTL_H)
	$(CPP) -c $(CPPFLAGS) $< -o $@

$(SHODIR)/%.sho : %.cc $(FTL_H)
	$(CPP) -c $(CPPSHFLAGS) $< -o $@

t%: $(ODIR)/t%.o lib/libFilter.a 
	$(LD) $(LDFLAGS) $< -o $@ $(LIBOPT)

lib: lib/libFilter.a

lib/libFilter.a: $(FTL_OBJ)
	$(AR) $@ $(FTL_OBJ) 

shared: lib/Filter.so

lib/Filter.so: $(FTL_SHO)
	\rm -rf lib/Filter.so
	$(LD) $(SHLDFLAGS) -o $(WDIR)/lib/Filter.so $(FTL_SHO)

Filter_dict.cc: $(WSH_H)
	rootcint -f Filter_dict.cc -c -p $(CINTOPT) $(WSH_H)  

mkfcomplex.cc: $(WSH_H)
	gcc -fPIC  -o shared/mkfcomplex.sho -c mkfcomplex.cc
	gcc -fPIC  -o obj/mkfcomplex.o -c mkfcomplex.cc

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
tar: archive/Filter.tar archive/macro.tar

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
	Filter_LinkDef.h Makefile INSTALL | \
        psnup -2 > archive/Filter.ps

