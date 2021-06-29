SETUP :
	/bin/tcsh ${HOME_CWB}/scripts/cwb_mkdir.csh .
	cp ${HOME_WAT}/tools/cwb/examples/S6A_R4_BKG_L1H1V1/input/* input/ 
	cp ${HOME_WAT}/tools/cwb/plugins/CWB_Plugin_HEN_SIM.C macro/ 
	cp ${HOME_WAT}/tools/cwb/plugins/CWB_Plugin_HEN_SIM_Config.C macro/ 

CONFIG :
	root -l -b macro/CreateSegmentHEN_S6A_L1H1V1.C

JOB :
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 1

CONDOR :
	root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_MACROS}/cwb_condor_create.C

SUBMIT :
	root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_MACROS}/cwb_condor_submit.C

MERGE :
	/bin/tcsh ${HOME_CWB}/scripts/cwb_merge.csh M1 

REPORT :
	/bin/tcsh ${HOME_CWB}/scripts/cwb_report.csh M1 create

CLEAN : 
	if [ -d macro ]; then rm macro/CWB_Plugin_HEN_SIM.C ; fi;
	if [ -d macro ]; then rm macro/CWB_Plugin_HEN_SIM_Config.C ; fi;
	if [ -d condor ]; then rm -rf condor ; fi;
	if [ -d input ]; then rm -rf input ; fi;
	if [ -d log ]; then rm -rf log ; fi;
	if [ -d output ]; then rm -rf output ; fi;
	if [ -d merge ]; then rm -rf merge ; fi;
	if [ -d report ]; then rm -rf report ; fi;
	if [ -d data ]; then rm -rf data ; fi;
	if [ -d tmp ]; then rm -rf tmp ; fi;
