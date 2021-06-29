setup :
	/bin/tcsh ${HOME_CWB}/scripts/cwb_mkdir.csh .
	if [ -d plugins ]; then rm plugins ; fi;
	ln -s ${HOME_CWB}/plugins plugins
	if [ ! -d frames ]; then mkdir frames ; fi;
	root -l -b macro/CreateFramesNOISE.C
	ls -1 ${PWD}/frames/L1H1V1-SimStrain-9311/H1*.gwf > input/H1.frames
	ls -1 ${PWD}/frames/L1H1V1-SimStrain-9311/L1*.gwf > input/L1.frames
	ls -1 ${PWD}/frames/L1H1V1-SimStrain-9311/V1*.gwf > input/V1.frames
	if [ -d input/inspiral.in ]; then rm input/inspiral.in ; fi;
	echo "931158200    931158600" > input/inspiral.in

cstrain :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data/*
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet2G.csh config/user_parameters.C CSTRAIN 1

coherence :
	if [ ! -d frames ]; then make setup ; fi;
	rm -f data/wave*
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet2G.csh data/cstrain_931158208_192_ADV_SIM_NSNS_L1H1V1_2G_ChirpletChains_60_job1.root COHERENCE config/user_parameters.C

super :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data/*
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet2G.csh config/user_parameters.C SUPERCLUSTER 1

like :
	if [ ! -d frames ]; then make setup ; fi;
	rm -f data/wave*
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet2G.csh data/supercluster_931158208_192_ADV_SIM_NSNS_L1H1V1_2G_ChirpletChains_60_job1.root LIKELIHOOD config/sparse_data_parameters.C

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
