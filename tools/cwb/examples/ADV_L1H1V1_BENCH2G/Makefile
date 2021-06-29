
ifndef MODE
MODE = SGQ9
endif

help : 
	@echo ''
	@echo 'make setup MODE=SGQ9_SSRA  : SGQ9 simulation with SSRA' 
	@echo 'make setup MODE=BBH_gMRA   : BBH  simulation with gMRA'
	@echo 'make setup MODE=BNS_pMRA   : BNS  simulation with pMRA'
	@echo 'make setup MODE=BKG_RSRA   :      background with RSRA'
	@echo ''
	@echo 'make cwb                   : start job=25 in interactive mode'
	@echo 'make ced                   : start job=25 in interactive mode'
	@echo 'make condor                : create batch condor files'
	@echo 'make batch                 : start batch condor files'
	@echo 'make merge                 : merge output root files'
	@echo 'make report                : create report'
	@echo 'make clean                 : delete files produced by the analysis'
	@echo 'make reset                 : remove all files created by setup and analysis'
	@echo ''

ifeq (${MODE},BKG)
CPLUGIN	= root -l -b -q macro/CWB_Plugin_SimNoise.C++
else
CPLUGIN  = root -l -b -q macro/CWB_Plugin_SimNoise_MDC_OTF.C++
endif

setup :	clean
	/bin/tcsh ${HOME_CWB}/scripts/cwb_mkdir.csh .
	if [ -d plugins ]; then rm plugins ; fi;
	ln -s ${HOME_CWB}/plugins plugins
	cd macro; ln -sf ../plugins/CWB_Plugin_SimNoise.C .; cd ..
	if [ -d input/period.in ]; then rm input/period.in ; fi;
	echo "1096827000      1096967000" > input/period.in
	cd config; ln -sf user_parameters_${MODE}.C  user_parameters.C; cd ..
	cd config; ln -sf user_pparameters_${MODE}.C user_pparameters.C; cd ..
	/bin/tcsh ${HOME_CWB}/scripts/cwb_mplugin.csh macro/CWB_Plugin_SimNoise_MDC_OTF.C \
                              plugins/CWB_Plugin_SimNoise.C plugins/CWB_Plugin_MDC_OTF.C
	${CPLUGIN}

cwb :	clean
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 25

ced : 	clean
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 25 0 true

condor : 
	/bin/tcsh ${HOME_CWB}/scripts/cwb_condor.csh create

bactch : 
	/bin/tcsh ${HOME_CWB}/scripts/cwb_condor.csh submit

merge :
	/bin/tcsh ${HOME_CWB}/scripts/cwb_merge.csh M1

report :
	/bin/tcsh ${HOME_CWB}/scripts/cwb_report.csh M1 create

clean : 
	if [ -d condor ]; then rm -f condor/*.dag*;rm -f condor/*.sub*; fi;
	if [ -d log ]; then rm -f log/* ; fi;
	if [ -d output ]; then rm -f output/* ; fi;
	if [ -d merge ]; then rm -f merge/* ; fi;
	if [ -d report/postprod ]; then rm -rf report/postprod/* ; fi;
	if [ -d report/dump ]; then rm -rf report/dump/* ; fi;
	if [ -d report/ced ]; then rm -rf report/ced/* ; fi;
	if [ -d data ]; then rm -rf data/* ; fi;
	if [ -d tmp ]; then rm -rf tmp/* ; fi;

reset :
	if [ -d plugins ]; then rm  plugins ; fi;
	if [ -d condor ]; then rm -rf condor ; fi;
	if [ -d log ]; then rm -rf log ; fi;
	if [ -d output ]; then rm -rf output ; fi;
	if [ -d merge ]; then rm -rf merge ; fi;
	if [ -d report ]; then rm -rf report ; fi;
	if [ -d data ]; then rm -rf data ; fi;
	if [ -d tmp ]; then rm -rf tmp ; fi;
	rm -f macro/CWB_Plugin_SimNoise_MDC_OTF.C
	rm -f macro/CWB_Plugin_SimNoise_MDC_OTF_C.d
	rm -f macro/CWB_Plugin_SimNoise_MDC_OTF_C.so
	rm -f macro/CWB_Plugin_SimNoise.C
	rm -f macro/CWB_Plugin_SimNoise_C.d
	rm -f macro/CWB_Plugin_SimNoise_C.so

