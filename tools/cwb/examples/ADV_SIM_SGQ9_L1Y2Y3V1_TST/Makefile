setup :
	/bin/tcsh ${HOME_CWB}/scripts/cwb_mkdir.csh batch .
	if [ -d plugins ]; then rm plugins ; fi;
	ln -s ${HOME_CWB}/plugins plugins
	if [ ! -d frames ]; then mkdir frames ; fi;
	root -l -b macro/CreateFramesNOISE.C
	root -l -b macro/CreateFramesMDC.C	 
	ls -1 ${PWD}/frames/L1Y2Y3V1-SimStrain-9311/Y2*.gwf > input/Y2.frames
	ls -1 ${PWD}/frames/L1Y2Y3V1-SimStrain-9311/Y3*.gwf > input/Y3.frames
	ls -1 ${PWD}/frames/L1Y2Y3V1-SimStrain-9311/L1*.gwf > input/L1.frames
	ls -1 ${PWD}/frames/L1Y2Y3V1-SimStrain-9311/V1*.gwf > input/V1.frames
	ls -1 ${PWD}/frames/L1Y2Y3V1-SGQ9-9311/*.gwf > input/MDC.frames
	if [ -d input/burst.in ]; then rm input/burst.in ; fi;
	echo "931158370    931158430" > input/burst.in

cwb :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf data/*
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 1 931158395 false

ced :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf data/*
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 1 931158395 true

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
