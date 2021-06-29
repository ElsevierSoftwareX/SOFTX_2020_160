WAT      = wat
TOOLS    = tools

ifndef ODIR
ODIR        = "."
endif

ifndef XIFO
XIFO        = 4
endif

all: all_wat all_tools html_footer html_header

all_wat:
ifdef WAT
	cd wat ; make XIFO=${XIFO}
endif
all_tools:
ifdef TOOLS
	cd tools ; make XIFO=${XIFO}
endif
html_footer:
	cd html ; make html_footer
html_header:
	cd html ; make html_header

sl_syspkg:
	make -f Makefile.sl_syspkg
	
deb_syspkg:
	make -f Makefile.deb_syspkg
	
clean: clean_wat clean_tools clean_html

clean_wat:
ifdef WAT
	cd wat ; make clean
endif
clean_tools:
ifdef TOOLS
	cd tools ; make clean
endif
clean_html:
	cd html ; make clean

install: winstall
uninstall: wuninstall

winstall:
	if [ ${ODIR} = "." ]; then 						\
	  cd tools ; make winstall ODIR=${ODIR}; 				\
	else 									\
	  if [ -d ${ODIR} ]; then 						\
	    echo "";								\
	    echo "Error: install directory already exist. Install aborted!!!"; 	\
	    echo "";								\
	  else 									\
	    cd tools ; make winstall ODIR=${ODIR}; 				\
	  fi; 									\
	fi

wuninstall:
	cd tools ; make wuninstall

ALL : all winstall
