
# ------------------------------------------------------------------------------
# signals
# ------------------------------------------------------------------------------

setup :
	/bin/tcsh ${HOME_CWB}/scripts/cwb_mkdir.csh .
	ln -sf ${HOME_CWB}/plugins plugins
	if [ ! -d frames ]; then mkdir frames ; fi;
	if [ ! -d frames/logs ]; then mkdir frames/logs ; fi;
	if [ -d input/frame.in ]; then rm input/frame.in ; fi;
	echo "1040738992        1040800000" > input/frame.in
	if [ -d input/burst.in ]; then rm input/burst.in ; fi;
	echo "1040738992        1040800000" > input/burst.in
	export CWB_UPARAMETERS_FILE="./config/frame_parameters.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_dump.csh sjob		
	/bin/tcsh ./edgesDQ.csh 

batch :
	rm -rf tmp/*
	rm -rf frames/*
	mkdir frames/logs
	export CWB_UPARAMETERS_FILE="./config/frame_parameters.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_condor.csh create

submit :
	export CWB_UPARAMETERS_FILE="./config/frame_parameters.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_condor.csh submit

mlog :
	export CWB_UPARAMETERS_FILE="./config/frame_parameters.C" ;	\
	/bin/tcsh ./mlog.csh 

list :
	export CWB_UPARAMETERS_FILE="./config/frame_parameters.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_dump.csh sjob		

dq :
	export CWB_UPARAMETERS_FILE="./config/frame_parameters.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_dump.csh dq		
	/bin/tcsh ./edgesDQ.csh 

frame :
	rm -rf tmp/*
	rm -rf frames/*
	mkdir -p frames/logs
	export CWB_UPARAMETERS_FILE="./config/frame_parameters.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 2 0 false 0 0 frames
	ls -1 ${PWD}/frames/L1/*.gwf > input/L1.lst
	ls -1 ${PWD}/frames/H1/*.gwf > input/H1.lst
	ls -1 ${PWD}/frames/V1/*.gwf > input/V1.lst

sim :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data/*
	export CWB_UPARAMETERS_FILE="./config/sim_parameters.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 2

sced :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data/*
	export CWB_UPARAMETERS_FILE="./config/sim_parameters.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 2 0 true

zprod :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data/*
	export CWB_UPARAMETERS_FILE="./config/zprod_parameters.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 2

zpced :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data/*
	export CWB_UPARAMETERS_FILE="./config/zprod_parameters.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 2 0 true

mzprod :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data/*
	export CWB_UPARAMETERS_FILE="./config/mzprod_parameters.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 2

mzpced :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data/*
	export CWB_UPARAMETERS_FILE="./config/mzprod_parameters.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 2 0 true

# ------------------------------------------------------------------------------
# signals+glitches
# ------------------------------------------------------------------------------

gsetup :
	/bin/tcsh ${HOME_CWB}/scripts/cwb_mkdir.csh 
	ln -sf ${HOME_CWB}/plugins plugins
	if [ ! -d frames ]; then mkdir frames ; fi;
	if [ ! -d frames/logs ]; then mkdir frames/logs ; fi;
	if [ -d input/frame.in ]; then rm input/frame.in ; fi;
	echo "1040738992        1040800000" > input/frame.in
	if [ -d input/burst.in ]; then rm input/burst.in ; fi;
	echo "1040738992        1040800000" > input/burst.in
	export CWB_UPARAMETERS_FILE="./config/gframe_parameters.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_dump.csh sjob		
	/bin/tcsh ./edgesDQ.csh 


gbatch :
	rm -rf tmp/*
	rm -rf frames/*
	mkdir frames/logs
	export CWB_UPARAMETERS_FILE="./config/gframe_parameters.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_condor.csh create

gsubmit :
	export CWB_UPARAMETERS_FILE="./config/gframe_parameters.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_condor.csh submit

gmlog :
	export CWB_UPARAMETERS_FILE="./config/gframe_parameters.C" ;	\
	/bin/tcsh ./gmlog.csh 

glist :
	export CWB_UPARAMETERS_FILE="./config/gframe_parameters.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_dump.csh sjob		

gdq :
	export CWB_UPARAMETERS_FILE="./config/gframe_parameters.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_dump.csh sjob		
	/bin/tcsh ./edgesDQ.csh 

glitch :
	root -l -b macro/CreateListOfGlitches.C

gframe :
	rm -rf tmp/*
	rm -rf frames/*
	mkdir -p frames/logs
	export CWB_UPARAMETERS_FILE="./config/gframe_parameters.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 2 0 false 0 0 frames
	ls -1 ${PWD}/frames/L1/*.gwf > input/L1.lst
	ls -1 ${PWD}/frames/H1/*.gwf > input/H1.lst 
	ls -1 ${PWD}/frames/V1/*.gwf > input/V1.lst 

gsim :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data/*
	export CWB_UPARAMETERS_FILE="./config/gsim_parameters.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 2

gsced :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data/*
	export CWB_UPARAMETERS_FILE="./config/gsim_parameters.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 2 0 true

gzprod :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data/*
	export CWB_UPARAMETERS_FILE="./config/gzprod_parameters.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 2

gzpced :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data/*
	export CWB_UPARAMETERS_FILE="./config/gzprod_parameters.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 2 0 true

gmzprod :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data/*
	export CWB_UPARAMETERS_FILE="./config/gmzprod_parameters.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 2

gmzpced :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data/*
	export CWB_UPARAMETERS_FILE="./config/gmzprod_parameters.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 2 0 true

clean : 
	if [ -d plugins ]; then rm  plugins ; fi;
	if [ -d frames ]; then rm -rf frames ; fi;
	if [ -d condor ]; then rm -rf condor ; fi;
	if [ -d log ]; then rm -rf log ; fi;
	if [ -d output ]; then rm -rf output ; fi;
	if [ -d merge ]; then rm -rf merge ; fi;
	if [ -d report ]; then rm -rf report ; fi;
	if [ -d data ]; then rm -rf data ; fi;
	if [ -d input ]; then rm -rf input ; fi;
	if [ -d tmp ]; then rm -rf tmp ; fi;
