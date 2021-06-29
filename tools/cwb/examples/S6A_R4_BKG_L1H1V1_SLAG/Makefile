setup :
	/bin/tcsh ${HOME_CWB}/scripts/cwb_mkdir.csh .

cwb :
	rm -rf tmp/*
	rm -rf data/*
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 8

ced :
	rm -rf tmp/*
	rm -rf data/*
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 8 0 true

clean : 
	if [ -d condor ]; then rm -rf condor ; fi;
	if [ -d log ]; then rm -rf log ; fi;
	if [ -d output ]; then rm -rf output ; fi;
	if [ -d merge ]; then rm -rf merge ; fi;
	if [ -d report ]; then rm -rf report ; fi;
	if [ -d data ]; then rm -rf data ; fi;
	if [ -d macro ]; then rm -rf macro ; fi;
	if [ -d tmp ]; then rm -rf tmp ; fi;
