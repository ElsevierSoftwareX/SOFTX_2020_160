setup :
	/bin/tcsh ${HOME_CWB}/scripts/cwb_mkdir.csh .
	if [ -d plugins ]; then rm plugins ; fi;
	ln -s ${HOME_CWB}/plugins plugins
	cp ${HOME_CWB}/plugins/CWB_Plugin_MDC_OTF.C macro/.
	if [ ! -d frames ]; then mkdir frames ; fi;
	root -l -b macro/CreateFramesNOISE.C
	ls -1 ${PWD}/frames/L1H1V1-SimStrain-9311/H1*.gwf > input/H1.frames
	ls -1 ${PWD}/frames/L1H1V1-SimStrain-9311/L1*.gwf > input/L1.frames
	ls -1 ${PWD}/frames/L1H1V1-SimStrain-9311/V1*.gwf > input/V1.frames
	if [ -d input/burst.in ]; then rm input/burst.in ; fi;
	echo "931158200    931158600" > input/inspiral.in

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
