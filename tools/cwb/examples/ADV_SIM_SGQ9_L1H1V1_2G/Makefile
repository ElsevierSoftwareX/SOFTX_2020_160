setup :
	/bin/tcsh ${HOME_CWB}/scripts/cwb_mkdir.csh  .
	if [ -d plugins ]; then rm plugins ; fi;
	ln -s ${HOME_CWB}/plugins plugins
	if [ ! -d frames ]; then mkdir frames ; fi;
	root -l -b macro/CreateFramesNOISE.C
	root -l -b macro/CreateFramesMDC.C	 
	ls -1 ${PWD}/frames/L1H1V1-SimStrain-9311/H1*.gwf > input/H1.frames
	ls -1 ${PWD}/frames/L1H1V1-SimStrain-9311/L1*.gwf > input/L1.frames
	ls -1 ${PWD}/frames/L1H1V1-SimStrain-9311/V1*.gwf > input/V1.frames
	ls -1 ${PWD}/frames/L1H1V1-SGQ9-9311/*.gwf > input/MDC.frames
	if [ -d input/burst.in ]; then rm input/burst.in ; fi;
	echo "931158360    931158440" > input/burst.in

cwb :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data/*
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 1 931158395 false

ced :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data/*
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 1 931158395 true
	mv data/ced_931158370_60_ADV_SIM_SGQ9_L1H1V1_2G_*_30_job1/ report/ced/

full :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data/*
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet2G.csh config/user_parameters.C FULL 1

init :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data/*
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet2G.csh config/user_parameters.C INIT 1

strain :
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet2G.csh data/init_931158378_44_ADV_SIM_SGQ9_L1H1V1_2G_job1.root STRAIN

cstrain :
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet2G.csh data/strain_931158378_44_ADV_SIM_SGQ9_L1H1V1_2G_job1.root CSTRAIN

coherence :
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet2G.csh data/cstrain_931158378_44_ADV_SIM_SGQ9_L1H1V1_2G_30_job1.root COHERENCE

supercluster :
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet2G.csh data/coherence_931158378_44_ADV_SIM_SGQ9_L1H1V1_2G_30_job1.root SUPERCLUSTER

likelihood :
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet2G.csh data/supercluster_931158378_44_ADV_SIM_SGQ9_L1H1V1_2G_30_job1.root LIKELIHOOD

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
