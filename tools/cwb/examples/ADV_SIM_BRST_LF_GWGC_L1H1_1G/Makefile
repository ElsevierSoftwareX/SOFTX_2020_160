SETUP :
	/bin/tcsh ${HOME_CWB}/scripts/cwb_mkdir.csh .
	ln -s ${HOME_CWB}/plugins plugins
	cp ${HOME_WAT}/tools/cwb/tutorials/RootSkyMap2Fits.C macro/ 
	cp ${HOME_WAT}/tools/cwb/tutorials/RootsSkyMap2Fits.C macro/ 
	cp ${HOME_WAT}/tools/cwb/tutorials/ReadSkyMapFromTree.C macro/ 
	cp ${HOME_WAT}/tools/gwat/tutorials/DrawFits.C macro/ 

JOB :
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 1

CED :
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 1 0 true

CONDOR :
	root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_MACROS}/cwb_condor_create.C

SUBMIT :
	root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_MACROS}/cwb_condor_submit.C

MERGE :
	/bin/tcsh ${HOME_CWB}/scripts/cwb_merge.csh M1 

REPORT :
	/bin/tcsh ${HOME_CWB}/scripts/cwb_report.csh M1 create

CLEAN : 
	if [ -d macro ]; then rm -f macro/RootSkyMap2Fits.C ; fi;
	if [ -d macro ]; then rm -f macro/RootsSkyMap2Fits.C ; fi;
	if [ -d macro ]; then rm -f macro/ReadSkyMapFromTree.C ; fi;
	if [ -d macro ]; then rm -f macro/DrawFits.C ; fi;
	if [ -d condor ]; then rm -rf condor ; fi;
	if [ -d plugins ]; then rm plugins ; fi;
	if [ -d log ]; then rm -rf log ; fi;
	if [ -d output ]; then rm -rf output ; fi;
	if [ -d merge ]; then rm -rf merge ; fi;
	if [ -d report ]; then rm -rf report ; fi;
	if [ -d data ]; then rm -rf data ; fi;
	if [ -d tmp ]; then rm -rf tmp ; fi;
