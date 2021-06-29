
ifndef FACTOR
FACTOR = 0
endif


setup :
	/bin/tcsh ${HOME_CWB}/scripts/cwb_mkdir.csh .
	if [ -d plugins ]; then rm plugins ; fi;
	ln -s ${HOME_CWB}/plugins plugins
	if [ ! -d frames ]; then mkdir frames ; fi;
	root -l -b macro/CreateFramesNOISE.C
	ls -1 ${PWD}/frames/L1H1H2H3V1-SimStrain-9311/H1*.gwf > input/H1.frames
	ls -1 ${PWD}/frames/L1H1H2H3V1-SimStrain-9311/H2*.gwf > input/H2.frames
	ls -1 ${PWD}/frames/L1H1H2H3V1-SimStrain-9311/H3*.gwf > input/H3.frames
	ls -1 ${PWD}/frames/L1H1H2H3V1-SimStrain-9311/L1*.gwf > input/L1.frames
	ls -1 ${PWD}/frames/L1H1H2H3V1-SimStrain-9311/V1*.gwf > input/V1.frames
	if [ -d input/burst.in ]; then rm input/burst.in ; fi;
	echo "931158370    931158430" > input/burst.in
	cp ${HOME_WAT}/tools/gwat/tutorials/Test_H2H3_45deg.C macro/

cwb_L1H1V1 :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data/*
	rm -rf report/dump/*
	export CWB_UPARAMETERS_FILE="config/user_parameters_L1H1V1.C" ;		\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 1 0 false ${FACTOR}

ced_L1H1V1 :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data/*
	rm -rf report/dump/*
	export CWB_UPARAMETERS_FILE="config/user_parameters_L1H1V1.C" ;		\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 1 0 true ${FACTOR}

cwb_L1H2H3V1 :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data/*	
	rm -rf report/dump/*
	export CWB_UPARAMETERS_FILE="config/user_parameters_L1H2H3V1.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 1 0 false ${FACTOR}

ced_L1H2H3V1 :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data/*
	rm -rf report/dump/*
	export CWB_UPARAMETERS_FILE="config/user_parameters_L1H2H3V1.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 1 0 true ${FACTOR}

cwb_H1 :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data/*	
	rm -rf report/dump/*
	export CWB_UPARAMETERS_FILE="config/user_parameters_H1.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 1 0 false ${FACTOR}

cwb_H2 :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data/*	
	rm -rf report/dump/*
	export CWB_UPARAMETERS_FILE="config/user_parameters_H2.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 1 0 false ${FACTOR}

cwb_H3 :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data/*	
	rm -rf report/dump/*
	export CWB_UPARAMETERS_FILE="config/user_parameters_H3.C" ;	\
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 1 0 false ${FACTOR}

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
