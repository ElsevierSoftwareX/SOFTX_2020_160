
help : 
	@echo ''
	@echo 'make setup                 : SGQ9 simulation with SSRA' 
	@echo ''
	@echo 'make cwb                   : start job=25 in interactive mode'
	@echo 'make ced                   : start job=25 in interactive mode'
	@echo 'make skymask               : create animated skymask plot in plot dir'
	@echo ''
	@echo 'make reset                 : remove all files created by setup and analysis'
	@echo ''

CPLUGIN  = root -l -b -q macro/CWB_Plugin_SimNoise_MDC_OTF.C++

setup :	clean
	/bin/tcsh ${HOME_CWB}/scripts/cwb_mkdir.csh .
	if [ -d plugins ]; then rm plugins ; fi;
	ln -s ${HOME_CWB}/plugins plugins
	cd macro; ln -sf ../plugins/CWB_Plugin_SimNoise.C .; cd ..
	if [ -d input/period.in ]; then rm input/period.in ; fi;
	echo "870772910      871215278" > input/period.in
	/bin/tcsh ${HOME_CWB}/scripts/cwb_mplugin.csh macro/CWB_Plugin_SimNoise_MDC_OTF.C \
                              plugins/CWB_Plugin_SimNoise.C plugins/CWB_Plugin_MDC_OTF.C  \
                              macro/CWB_Plugin_RingSkyMask.C
	${CPLUGIN}

cwb :	clean
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 25

ced : 	clean
	/bin/tcsh ${HOME_CWB}/scripts/cwb_inet.csh 25 0 true

skymask : 
	root -l -q macro/RingSkyMask.C

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
	if [ -d plot ]; then rm -rf plot/* ; fi;

reset :
	if [ -d plugins ]; then rm  plugins ; fi;
	if [ -d condor ]; then rm -rf condor ; fi;
	if [ -d log ]; then rm -rf log ; fi;
	if [ -d output ]; then rm -rf output ; fi;
	if [ -d merge ]; then rm -rf merge ; fi;
	if [ -d report ]; then rm -rf report ; fi;
	if [ -d data ]; then rm -rf data ; fi;
	if [ -d tmp ]; then rm -rf tmp ; fi;
	if [ -d plot ]; then rm -rf plot/* ; fi;
	rm -f macro/CWB_Plugin_SimNoise_MDC_OTF.C
	rm -f macro/CWB_Plugin_SimNoise_MDC_OTF_C.d
	rm -f macro/CWB_Plugin_SimNoise_MDC_OTF_C.so
	rm -f macro/CWB_Plugin_SimNoise.C
	rm -f macro/CWB_Plugin_SimNoise_C.d
	rm -f macro/CWB_Plugin_SimNoise_C.so

