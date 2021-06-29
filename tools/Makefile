CWBSTFT     = tools/stft
CWBTOOLBOX  = tools/toolbox
CWBHISTORY  = tools/history
CWBGWAT     = tools/gwat
CWBBICO     = tools/bico
CWBFILTER   = tools/filter
CWBFRAME    = tools/frame
CWBCWB      = tools/cwb
CWBWAVEGRAPH= tools/wavegraph
CWBFRDISPLAY= tools/frdisplay
CWBWDMUTILS = tools/wdmutils

ifndef ODIR
ODIR        = "."
endif

ifndef XIFO
XIFO        = 4
endif

all: all_stft all_toolbox all_history all_gwat all_bico all_cwb all_wavegraph all_filter all_frame all_cwb all_frdisplay all_ebbh all_wdmutils
	chmod +x *.sh;
	chmod +x *.csh;

all_stft:
ifdef CWBSTFT
	cd stft ; make XIFO=${XIFO}
endif
all_toolbox:
ifdef CWBTOOLBOX
	cd toolbox ; make XIFO=${XIFO}
endif
all_history:
ifdef CWBHISTORY
	cd history ; make XIFO=${XIFO}
endif
all_gwat:
ifdef CWBGWAT
	cd gwat ; make XIFO=${XIFO}
endif
all_bico:
ifdef CWBBICO
	cd bico ; make  XIFO=${XIFO}
endif
all_cwb:
ifdef CWBCWB 
	cd cwb ; make XIFO=${XIFO}
endif
all_wavegraph:
ifdef CWBWAVEGRAPH 
	cd wavegraph ; make XIFO=${XIFO}
endif
all_filter:
ifdef CWBFILTER
	cd filter ; make  XIFO=${XIFO}
endif
all_frame:
ifdef CWBFRAME 
	cd frame ; make XIFO=${XIFO}
endif
all_frdisplay:
ifdef CWBFRDISPLAY 
	cd frdisplay ; make XIFO=${XIFO}
endif
all_wdmutils:
ifdef CWBWDMUTILS 
	cd wdmutils ; make XIFO=${XIFO}
endif
all_ebbh:
ifdef _USE_EBBH 
	cd eBBH ; make 
endif
sh: 
	root -l -b cwb/macros/macros/csh2sh.C;

clean: clean_stft clean_toolbox clean_history clean_gwat clean_bico clean_cwb clean_wavegraph clean_filter clean_frame clean_frdisplay clean_ebbh clean_wdmutils uninstall

clean_stft:
ifdef CWBSTFT
	cd stft ; make clean
endif
clean_toolbox:
ifdef CWBTOOLBOX
	cd toolbox ; make clean
endif
clean_history:
ifdef CWBHISTORY
	cd history ; make clean
endif
clean_gwat:
ifdef CWBGWAT
	cd gwat ; make clean
endif
clean_bico:
ifdef CWBBICO
	cd bico ; make clean
endif
clean_cwb:
ifdef CWBCWB
	cd cwb ; make clean
endif
clean_wavegraph:
ifdef CWBWAVEGRAPH
	cd wavegraph ; make clean
endif
clean_filter:
ifdef CWBFILTER
	cd filter ; make clean
endif
clean_frame:
ifdef CWBFRAME
	cd frame ; make clean
endif
clean_frdisplay:
ifdef CWBFRDISPLAY
	cd frdisplay ; make clean
endif
clean_ebbh:
ifdef _USE_EBBH
	cd eBBH ; make clean
endif
clean_wdmutils:
ifdef CWBWDMUTILS
	cd wdmutils ; make clean
endif

install: winstall
uninstall: wuninstall

wuninstall:
	if [ ${ODIR} = "." ]; then 							\
	  make remove  WAT_INSTALL=$(HOME_WAT_INSTALL);					\
	fi

winstall:
	cd ${HOME_WAT};									\
	if [ ${ODIR} = "." ]; then 							\
  	  echo ""; 									\
  	  echo "Install libraries and includes locally"; 				\
	  echo ""; 									\
	  make -f tools/Makefile copy  WAT_INSTALL=$(HOME_WAT_INSTALL);			\
	else 										\
	  GIT_STATUS=$$(git -C ${HOME_WAT} status 2> /dev/null);			\
	  if [[ ! $${GIT_STATUS} =~ "working tree clean" ]]; then			\
	    echo ""; 									\
	    echo "Before installing the local git repository must be clean"; 		\
	    echo ""; 									\
	  else										\
	    if [ -d ${ODIR} ]; then 							\
	      echo ""; 									\
	      echo "Error: install directory already exist. Install aborted!!!"; 	\
	      echo ""; 									\
	    else 									\
	      git clone --single-branch --depth 1 file:///$(HOME_WAT) ${ODIR}; 		\
	      make -f tools/Makefile copy  WAT_INSTALL=$$(echo ${HOME_WAT_INSTALL} | sed 's|'${HOME_WAT}'|'${ODIR}'|g'); \
	      echo ""; 									\
	      echo "New cloned directory is: ${ODIR}"; 					\
	      echo ""; 									\
	    fi; 									\
	  fi;										\
	fi

remove:
	if [ -d $(WAT_INSTALL)/lib ];     then rm -f $(WAT_INSTALL)/lib/*     ; fi;
	if [ -d $(WAT_INSTALL)/inc ];     then rm -f $(WAT_INSTALL)/inc/*     ; fi;
	if [ -d $(WAT_INSTALL)/macros ];  then rm -f $(WAT_INSTALL)/macros/*  ; fi;
	if [ -d $(WAT_INSTALL)/scripts ]; then rm -f $(WAT_INSTALL)/scripts/* ; fi;

copy:
	cd ${HOME_WAT};								\
	if [ ! -d $(WAT_INSTALL)/lib ]; then mkdir -p $(WAT_INSTALL)/lib; fi;	\
	if [ ! -d $(WAT_INSTALL)/inc ]; then mkdir -p $(WAT_INSTALL)/inc; fi;	\
	cp $(HOME_WAT)/wat/lib/*.so 		$(WAT_INSTALL)/lib/;		\
	cp $(HOME_WAT)/tools/history/lib/*.so	$(WAT_INSTALL)/lib/;		\
	cp $(HOME_WAT)/tools/toolbox/lib/*.so 	$(WAT_INSTALL)/lib/;		\
	cp $(HOME_WAT)/tools/stft/lib/*.so 	$(WAT_INSTALL)/lib/;		\
	cp $(HOME_WAT)/tools/bico/lib/*.so 	$(WAT_INSTALL)/lib/;		\
	cp $(HOME_WAT)/tools/cwb/lib/*.so 	$(WAT_INSTALL)/lib/;		\
	cp $(HOME_WAT)/tools/wavegraph/lib/*.so	$(WAT_INSTALL)/lib/;		\
	cp $(HOME_WAT)/tools/filter/lib/*.so 	$(WAT_INSTALL)/lib/;		\
	cp $(HOME_WAT)/tools/frame/lib/*.so	$(WAT_INSTALL)/lib/;		\
	cp $(HOME_WAT)/tools/gwat/lib/*.so	$(WAT_INSTALL)/lib/;		\
	cp $(HOME_WAT)/wat/*.hh 		$(WAT_INSTALL)/inc/;		\
	cp $(HOME_WAT)/tools/history/*.hh	$(WAT_INSTALL)/inc/;		\
	cp $(HOME_WAT)/tools/toolbox/*.hh 	$(WAT_INSTALL)/inc/;		\
	cp $(HOME_WAT)/tools/stft/*.hh 		$(WAT_INSTALL)/inc/;		\
	cp $(HOME_WAT)/tools/bico/*.hh 		$(WAT_INSTALL)/inc/;		\
	cp $(HOME_WAT)/tools/cwb/*.hh 		$(WAT_INSTALL)/inc/;		\
	cp $(HOME_WAT)/tools/cwb/plugins/*.h 	$(WAT_INSTALL)/inc/;		\
	cp $(HOME_WAT)/tools/wavegraph/*.hh 	$(WAT_INSTALL)/inc/;		\
	cp $(HOME_WAT)/tools/filter/*.hh	$(WAT_INSTALL)/inc/;		\
	cp $(HOME_WAT)/tools/frame/*.hh		$(WAT_INSTALL)/inc/;		\
	cp $(HOME_WAT)/tools/gwat/*.hh		$(WAT_INSTALL)/inc/
ifdef _USE_ROOT6 									
	cd ${HOME_WAT};								\
	cp $(HOME_WAT)/wat/*.pcm 		$(WAT_INSTALL)/lib/;		\
	cp $(HOME_WAT)/tools/*/*.pcm		$(WAT_INSTALL)/lib/
endif 
ifdef _USE_EBBH	
	cd ${HOME_WAT};								\
	cp $(HOME_WAT)/tools/eBBH/lib/*.so	$(WAT_INSTALL)/lib/
endif 
ifdef _USE_EBBH	
	cd ${HOME_WAT};								\
	cp $(HOME_WAT)/tools/eBBH/*.hh		$(WAT_INSTALL)/inc/
endif
