setup :
	/bin/tcsh ${HOME_CWB}/scripts/cwb_mkdir.csh  .
	if [ -d plugins ]; then rm plugins ; fi;
	ln -s ${HOME_CWB}/plugins plugins
	mkdir data1 data2
	if [ ! -d frames ]; then mkdir frames ; fi;
	root -l -b macro/CreateFramesNOISE.C
	ls -1 ${PWD}/frames/L1H1V1-SimStrain-9311/H1*.gwf > input/H1.frames
	ls -1 ${PWD}/frames/L1H1V1-SimStrain-9311/L1*.gwf > input/L1.frames
	ls -1 ${PWD}/frames/L1H1V1-SimStrain-9311/V1*.gwf > input/V1.frames
	if [ -d input/burst.in ]; then rm input/burst.in ; fi;
	echo "931158370    931158430" > input/burst.in

cwb :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data/*
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet2G.csh config/user_parameters.C 0 931158395 data

cwb1 :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data1/*
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet2G.csh config/user_parameters.C 1 931158395 data1
	/bin/tcsh ${HOME_CWB}/scripts/cwb_merge_logmdc.csh 

cwb2 :
	if [ ! ${IF} ]; then echo -ne "\nBad input : make cwb2 IF=XXX\n" ; false ; fi; 
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data2/*
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet2G.csh ${IF} 2 931158395 data2

clean : 
	if [ -d plugins ]; then rm  plugins ; fi;
	if [ -d frames ]; then rm -rf frames ; fi;
	if [ -d condor ]; then rm -rf condor ; fi;
	if [ -d log ]; then rm -rf log ; fi;
	if [ -d output ]; then rm -rf output ; fi;
	if [ -d merge ]; then rm -rf merge ; fi;
	if [ -d report ]; then rm -rf report ; fi;
	if [ -d data ]; then rm -rf data ; fi;
	if [ -d data1 ]; then rm -rf data1 ; fi;
	if [ -d data2 ]; then rm -rf data2 ; fi;
	if [ -d input ]; then rm -rf input ; fi;
	if [ -d tmp ]; then rm -rf tmp ; fi;
