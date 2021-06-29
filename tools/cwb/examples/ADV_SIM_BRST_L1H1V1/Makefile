setup :
	/bin/tcsh ${HOME_CWB}/scripts/cwb_mkdir.csh .
	if [ -d plugins ]; then rm plugins ; fi;
	ln -s ${HOME_CWB}/plugins plugins
	if [ -d input/burst.in ]; then rm input/burst.in ; fi;
	echo "93100000        93106000" > input/burst.in

condor :
	/bin/tcsh ${HOME_CWB}/scripts/cwb_condor.csh create 
	/bin/tcsh ${HOME_CWB}/scripts/cwb_condor.csh submit 
	
cwb :
	rm -rf tmp/*
	rm -rf data/*
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 1 

ced :
	if [ ! -d frames ]; then make setup ; fi;
	rm -rf tmp/*
	rm -rf data/*
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 1 0 true

clean : 
	if [ -d plugins ]; then rm  plugins ; fi;
	if [ -d condor ]; then rm -rf condor ; fi;
	if [ -d log ]; then rm -rf log ; fi;
	if [ -d output ]; then rm -rf output ; fi;
	if [ -d merge ]; then rm -rf merge ; fi;
	if [ -d report ]; then rm -rf report ; fi;
	if [ -d data ]; then rm -rf data ; fi;
	if [ -d input ]; then rm -rf input ; fi;
	if [ -d tmp ]; then rm -rf tmp ; fi;
