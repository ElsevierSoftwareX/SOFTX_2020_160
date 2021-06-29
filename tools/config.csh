#!/bin/tcsh

if (! $?SITE_CLUSTER ) then
  echo "cwb setup error ..."
  echo "try 'source tools/site_config.sh' (where SITE_CLUSTER=ATLAS/CIT/CNAF/CASCINA/GATECH/USER_DEFINED)"
  exit
endif

setenv ICC_INIT_SCRIPT "/opt/intel/2018/intel.csh"

if ("$SITE_CLUSTER" == ATLAS) then

  if (! $?PYTHONPATH ) then
    setenv PYTHONPATH ""
  endif
  if (! $?DYLD_LIBRARY_PATH ) then
    setenv DYLD_LIBRARY_PATH ""
  endif
  if ("$DYLD_LIBRARY_PATH" !~ *${PYTHONPATH}*) then
    setenv DYLD_LIBRARY_PATH ${PYTHONPATH}:${DYLD_LIBRARY_PATH}
  endif
  if ("$LD_LIBRARY_PATH" !~ "*${PYTHONPATH}*") then
    setenv LD_LIBRARY_PATH ${PYTHONPATH}:${LD_LIBRARY_PATH}
  endif
  # must be declared after DYLD_LIBRARY_PATH & LD_LIBRARY_PATH
  if ("$PYTHONPATH" !~ "*${ROOTSYS}/lib:*") then
    setenv PYTHONPATH ${ROOTSYS}/lib:${PYTHONPATH}
  endif
  if ("$PYTHONPATH" !~ "${HOME_LAL}/lib/python2.7/site-packages:") then
    setenv PYTHONPATH ${HOME_LAL}/lib/python2.7/site-packages:${PYTHONPATH}
  endif 
  if ("$LD_LIBRARY_PATH" !~ "*/usr/lib/*") then
    setenv LD_LIBRARY_PATH /usr/lib/:${LD_LIBRARY_PATH}
  endif

  setenv NODE_DATA_DIR       "/atlas/user/X_HOST_NAME/X_USER_NAME"  	# user node data name
#  setenv NODE_DATA_DIR       "/local/user/X_USER_NAME"                  # user node data name
#  setenv CONDOR_LOG_DIR      "/local/user/X_USER_NAME"                   # condor log
  setenv CONDOR_LOG_DIR      "/atlas/user/X_HOST_NAME/X_USER_NAME/condor" # condor log
  setenv WWW_PUBLIC_DIR      "X_HOME/WWW/LSC/reports"                   # public dir

  setenv ICC_INIT_SCRIPT 	"/opt/intel/2018/intel.csh"

endif

if ("$SITE_CLUSTER" == CIT) then

  if (! $?PYTHONPATH ) then
    setenv PYTHONPATH ""
  endif
  if (! $?DYLD_LIBRARY_PATH ) then
    setenv DYLD_LIBRARY_PATH ""
  endif
  if ("$DYLD_LIBRARY_PATH" !~ *${PYTHONPATH}*) then
    setenv DYLD_LIBRARY_PATH ${PYTHONPATH}:${DYLD_LIBRARY_PATH}
  endif
  if ("$LD_LIBRARY_PATH" !~ "*${PYTHONPATH}*") then
    setenv LD_LIBRARY_PATH ${PYTHONPATH}:${LD_LIBRARY_PATH}
  endif
  # must be declared after DYLD_LIBRARY_PATH & LD_LIBRARY_PATH
  if ("$PYTHONPATH" !~ "*${ROOTSYS}/lib:*") then
    setenv PYTHONPATH ${ROOTSYS}/lib:${PYTHONPATH}
  endif
  if ("$PYTHONPATH" !~ "${HOME_LAL}/lib64/python2.7/site-packages:") then
    setenv PYTHONPATH ${HOME_LAL}/lib64/python2.7/site-packages:${PYTHONPATH}
  endif 
  if ("$LD_LIBRARY_PATH" !~ "*/usr/lib/:/usr/lib64/*") then
    setenv LD_LIBRARY_PATH /usr/lib/:/usr/lib64/:${LD_LIBRARY_PATH}
  endif

#  setenv NODE_DATA_DIR       "/data/X_HOST_NAME/X_USER_NAME"  # user node data name
  setenv NODE_DATA_DIR       "/local/X_USER_NAME"	      # user node data name
  setenv CONDOR_LOG_DIR      "/usr1/X_USER_NAME"              # condor log
  setenv WWW_PUBLIC_DIR      "X_HOME/public_html/reports"     # public dir

  setenv ICC_INIT_SCRIPT     "/ldcg/intel/2020u2/compilers_and_libraries_2020.2.254/linux/bin/iccvars.sh"

endif

if ("$SITE_CLUSTER" == CNAF) then

  if ("$LD_LIBRARY_PATH" !~ "*/usr/lib/:/usr/lib64/*") then
    setenv LD_LIBRARY_PATH /usr/lib/:/usr/lib64/:${LD_LIBRARY_PATH}
  endif

  setenv NODE_DATA_DIR       "/home/VIRGO/X_USER_NAME"	      # user node data name
  setenv CONDOR_LOG_DIR      ""
  setenv WWW_PUBLIC_DIR      ""

endif

if ("$SITE_CLUSTER" == CASCINA) then

  if ("$LD_LIBRARY_PATH" !~ "*/usr/lib/:/usr/lib64/*") then
    setenv LD_LIBRARY_PATH /usr/lib/:/usr/lib64/:${LD_LIBRARY_PATH}
  endif

  setenv NODE_DATA_DIR       ""
  setenv CONDOR_LOG_DIR      ""
  setenv WWW_PUBLIC_DIR      ""

endif

if ("$SITE_CLUSTER" == GATECH) then

  if ("$LD_LIBRARY_PATH" !~ "*${METAIOLIB}*") then
    setenv LD_LIBRARY_PATH ${METAIOLIB}:${LD_LIBRARY_PATH}
  endif

  setenv NODE_DATA_DIR       "X_HOME"
  setenv CONDOR_LOG_DIR      "condor"
  setenv WWW_PUBLIC_DIR      ""

endif

if ("$SITE_CLUSTER" != CIT && "$SITE_CLUSTER" != ATLAS && "$SITE_CLUSTER" != CASCINA && "$SITE_CLUSTER" != CNAF && "$SITE_CLUSTER" != GATECH) then

  setenv NODE_DATA_DIR       ""
  setenv CONDOR_LOG_DIR      ""
  setenv WWW_PUBLIC_DIR      ""

endif

# define the batch system
if (! $?BATCH_SYSTEM ) then
  setenv BATCH_SYSTEM  ""
endif
if ("$BATCH_SYSTEM" != PEGASUS && "$BATCH_SYSTEM" != LSF && "$BATCH_SYSTEM" != OSG && "$BATCH_SYSTEM" != "") then
  setenv BATCH_SYSTEM  ""      # enable CONDOR
endif
if ("$BATCH_SYSTEM" == CONDOR || "$BATCH_SYSTEM" == "") then
  unsetenv _USE_PEGASUS
  unsetenv _USE_LSF
  unsetenv _USE_OSG
endif
if ("$BATCH_SYSTEM" == PEGASUS) then
  if (! $?CWB_PEGASUS_WATENV ) then
    echo ""
    echo "env CWB_PEGASUS_WATENV is not defined !"
    exit
  endif
  if (! $?CWB_PEGASUS_SITE ) then
    echo ""
    echo "env CWB_PEGASUS_SITE is not defined !"
    exit
  endif
  setenv   _USE_PEGASUS 1      # enable PEGASUS workflow (development version)
  unsetenv _USE_LSF
  unsetenv _USE_OSG
endif
if ("$BATCH_SYSTEM" == LSF) then
  if (! $?LSF_QUEUE ) then
    echo ""
    echo "env LSF_QUEUE is not defined !"
    exit
  endif
  unsetenv _USE_PEGASUS
  setenv   _USE_LSF 1          # enable LSF (development version)
  unsetenv _USE_OSG
endif
if ("$BATCH_SYSTEM" == OSG) then
  unsetenv _USE_PEGASUS
  unsetenv _USE_LSF
  setenv   _USE_OSG 1          # enable OSG (development version)
endif

if (! $?_USE_ICC ) then
   echo -n ""
else
  if ( ! -e $ICC_INIT_SCRIPT ) then
    echo ""
    echo -n $ICC_INIT_SCRIPT
    echo " not exist!"
    echo "intel compiler can not be used, disable ICC in watenv -> unsetenv _USE_ICC"
    exit
  else
    if ("$LD_LIBRARY_PATH" !~ "intel64") then
      if ("$SITE_CLUSTER" == ATLAS) then
        source $ICC_INIT_SCRIPT
      else 
        source $ICC_INIT_SCRIPT -arc intel64
      endif 
    endif 
#    if ( $? != 0) then
#      echo ""
#      echo "source" $ICC_INIT_SCRIPT "intel64 error !!!"
#      echo "intel compiler can not be used, disable ICC in watenv -> unsetenv _USE_ICC"
#      echo ""
#      exit
#    endif
  endif
endif

if (! $?CWB_ANALYSIS ) then
  echo ""
  echo "CWB_ANALYSIS env not defined in watenv setup"
  echo "add in watenv setup the following statement"
  echo ""
  echo "setenv CWB_ANALYSIS        '1G'       # 1G  analysis "
  echo "or"
  echo "setenv CWB_ANALYSIS        '2G'       # 2G  analysis "
  echo ""
  echo "use 'cwb_setpipe [1G/1g/1,2G/2g/2]' command to switch analysis type [1G/2G]"
  exit
endif
if (( ${CWB_ANALYSIS} != '1G' ) && ( ${CWB_ANALYSIS} != '2G' )) then
  echo -n "CWB_ANALYSIS="
  echo ${CWB_ANALYSIS}
  echo "CWB_ANALYSIS env must be 1G or 2G"
  echo "check watenv setup"
  exit
endif

setenv CWB_TOOLBOX           "${HOME_WAT}/tools/toolbox"
setenv CWB_SKYPLOT           "${HOME_WAT}/tools/skyplot"
setenv CWB_GWAT              "${HOME_WAT}/tools/gwat"
setenv CWB_EBBH              "${HOME_WAT}/tools/eBBH"
setenv CWB_HISTORY           "${HOME_WAT}/tools/history"
setenv CWB_STFT              "${HOME_WAT}/tools/stft"
setenv CWB_BICO              "${HOME_WAT}/tools/bico"
setenv CWB_FILTER            "${HOME_WAT}/tools/filter"
setenv CWB_FRAME             "${HOME_WAT}/tools/frame"
setenv CWB_WAVESCAN          "${HOME_WAT}/tools/wavescan"
setenv CWB_ONLINE            "${HOME_WAT}/tools/online"

if ("$PATH" !~ "*${ROOTSYS}/bin*") then
  setenv PATH "${ROOTSYS}/bin:${PATH}"
endif
if ( "$OSTYPE" =~ "linux*" ) then
  if ("$PATH" !~ "*${HOME_FRLIB}/Linux*") then
    setenv PATH "${HOME_FRLIB}/Linux:${PATH}"
  endif
endif
if ( "$OSTYPE" =~ "darwin*" ) then
  if ("$PATH" !~ "*${HOME_FRLIB}/Darwin*") then
    setenv PATH "${HOME_FRLIB}/Darwin:${PATH}"
  endif
endif
if ("$LD_LIBRARY_PATH" !~ "*${ROOTSYS}/lib:*") then
  setenv LD_LIBRARY_PATH "${ROOTSYS}/lib:${LD_LIBRARY_PATH}"
endif
if ( $?_USE_ROOT6 ) then
  if ("$LD_LIBRARY_PATH" !~ "*${HOME_WAT_INSTALL}/lib:*") then
    setenv LD_LIBRARY_PATH "${HOME_WAT_INSTALL}/lib:${LD_LIBRARY_PATH}"
  endif
endif
if ( $?_USE_LAL ) then
  if ("$LD_LIBRARY_PATH" !~ "*${HOME_LAL}/lib:*") then
    setenv LD_LIBRARY_PATH "${HOME_LAL}/lib:${LD_LIBRARY_PATH}"
  endif
endif

# setenv ROOT_INCLUDE_PATH: used for include path by macros
setenv ROOT_INCLUDE_PATH "$CWB_CONFIG"

if (! $?MANPATH ) then
  setenv MANPATH ""
endif
if ("$MANPATH" !~ "*${HOME_WAT}/man*") then
  setenv MANPATH "${HOME_WAT}/man:${MANPATH}"
endif

alias  root                  "${ROOTSYS}/bin/root -l" 
alias  broot		     "${ROOTSYS}/bin/root -b -l"

setenv HOME_CWB              "${HOME_WAT}/tools/cwb"  
setenv CWB_SCRIPTS           "${HOME_CWB}/scripts"
setenv CWB_DIR_NAME          "${HOME_CWB}/macros"	# obsolete (only for back compatibility)
setenv CWB_MACROS            "${HOME_CWB}/macros"

setenv HOME_WAVEBICO         "${HOME_WAT}/tools/wavebico"  

setenv HOME_WAVEMDC          "${HOME_WAT}/tools/wavemdc"  
setenv WMDC_CONFIG           "./wmdc_config.C"  
alias  wmdc_condor           "${HOME_WAVEMDC}/wmdc_condor.csh" 		#_SKIP_csh2sh_

setenv HOME_FRDISPLAY        "${HOME_WAT}/tools/frdisplay"

setenv CWB_ROOTLOGON_FILE    "${CWB_MACROS}/cwb_rootlogon.C"

if ( "$CWB_ANALYSIS" == "1G" ) then
  setenv CWB_PARAMETERS_FILE   "${CWB_MACROS}/cwb1G_parameters.C"
endif
if ( "$CWB_ANALYSIS" == "2G" ) then
  setenv CWB_PARAMETERS_FILE   "${CWB_MACROS}/cwb2G_parameters.C"
endif
if (! $?CWB_UPARAMETERS_FILE ) then
  setenv CWB_UPARAMETERS_FILE  "./config/user_parameters.C"
endif
if (! $?CWB_MPARAMETERS_FILE ) then
  setenv CWB_MPARAMETERS_FILE  "./merge/cwb_user_parameters.C"
endif
if (! $?CWB_EMPARAMETERS_FILE ) then
  setenv CWB_EMPARAMETERS_FILE  ""
endif
setenv CWB_EPARAMETERS_FILE  "${CWB_MACROS}/cwb_eparameters.C"
setenv CWB_PPARAMETERS_FILE  "${CWB_MACROS}/cwb_pparameters.C"
if (! $?CWB_UPPARAMETERS_FILE ) then
  setenv CWB_UPPARAMETERS_FILE "./config/user_pparameters.C"
endif
setenv CWB_EPPARAMETERS_FILE "${CWB_MACROS}/cwb_epparameters.C"

setenv CWB_NETS_FILE         "${CWB_SCRIPTS}/cwb_net.sh"

if (! $?_USE_OSG ) then
  echo -n ""
else
  setenv CWB_NETS_FILE       "${CWB_SCRIPTS}/cwb_net_osg.sh"
endif
if (! $?_USE_LSF ) then
  echo -n ""
else
  setenv CWB_NETS_FILE       "${CWB_SCRIPTS}/cwb_net_lsf.sh"
endif
if (! $?_USE_PEGASUS ) then
  echo -n ""
else
  setenv CWB_NETS_FILE       "${CWB_SCRIPTS}/cwb_net_pegasus.sh"
endif

setenv CWB_NETC_FILE         "${CWB_MACROS}/cwb_net.C"

setenv CWB_HTML_INDEX        "${CWB_MACROS}/html_templates/html_index_template.txt"
setenv CWB_HTML_HEADER       "${CWB_MACROS}/html_templates/html_header_template.txt"
setenv CWB_HTML_BODY_PROD    "${CWB_MACROS}/html_templates/html_body_prod_template.txt"
setenv CWB_HTML_BODY_SIM     "${CWB_MACROS}/html_templates/html_body_sim_template.txt"

setenv CWB_LAG_NUMBER        -1
setenv CWB_SLAG_NUMBER       -1

alias cwb_setpipe            "source ${CWB_SCRIPTS}/cwb_setpipe.csh"
alias cwb_setpars            "source ${CWB_SCRIPTS}/cwb_setpars.csh" 
alias cwb_setppars           "source ${CWB_SCRIPTS}/cwb_setppars.csh" 
alias cwb_clonedir           "${CWB_SCRIPTS}/cwb_clonedir.csh" 		#_SKIP_csh2sh_

alias cwb_mkdir              "${CWB_SCRIPTS}/cwb_mkdir.csh" 		#_SKIP_csh2sh_
alias cwb_online             "${CWB_SCRIPTS}/cwb_online.csh" 		#_SKIP_csh2sh_
alias cwb_json2dat           "${CWB_SCRIPTS}/cwb_json2dat.csh" 		#_SKIP_csh2sh_
alias cwb_petools            "${CWB_SCRIPTS}/cwb_petools.csh" 		#_SKIP_csh2sh_
alias cwb_pereport           "${CWB_SCRIPTS}/cwb_pereport.csh" 		#_SKIP_csh2sh_

alias cwb_tune_slag          "${CWB_SCRIPTS}/cwb_tune_slag.csh" 	#_SKIP_csh2sh_
alias cwb_get_rejslag        "${CWB_SCRIPTS}/cwb_get_rejslag.csh" 	#_SKIP_csh2sh_

alias cwb_sub_lsf            "${CWB_SCRIPTS}/cwb_sub_lsf.csh" 		#_SKIP_csh2sh_

alias cwb_inet               "${CWB_SCRIPTS}/cwb_inet.csh" 		#_SKIP_csh2sh_
alias cwb_xnet               "${CWB_SCRIPTS}/cwb_xnet.csh" 		#_SKIP_csh2sh_
alias cwb_jnet               "${CWB_SCRIPTS}/cwb_jnet.csh" 		#_SKIP_csh2sh_
alias cwb_eced               "${CWB_SCRIPTS}/cwb_eced.csh" 		#_SKIP_csh2sh_
alias cwb_eced2G             "${CWB_SCRIPTS}/cwb_eced2G.csh" 		#_SKIP_csh2sh_
alias cwb_inet2G             "${CWB_SCRIPTS}/cwb_inet2G.csh" 		#_SKIP_csh2sh_
alias cwb_inet2G1            "${CWB_SCRIPTS}/cwb_inet2G1.csh" 		#_SKIP_csh2sh_
alias cwb_inet2G2            "${CWB_SCRIPTS}/cwb_inet2G2.csh" 		#_SKIP_csh2sh_
alias cwb_mplugin            "${CWB_SCRIPTS}/cwb_mplugin.csh" 		#_SKIP_csh2sh_
alias cwb_gwosc              "make -f ${HOME_CWB}/gwosc/Makefile.gwosc" 
alias cwb_compile            "${CWB_SCRIPTS}/cwb_compile.csh" 		#_SKIP_csh2sh_

alias cwb_merge              "${CWB_SCRIPTS}/cwb_merge.csh" 		#_SKIP_csh2sh_
alias cwb_xtalk              "${CWB_SCRIPTS}/cwb_xtalk.csh" 		#_SKIP_csh2sh_

alias cwb_publish            "${CWB_SCRIPTS}/cwb_publish.csh" 		#_SKIP_csh2sh_
alias cwb_publish_ced        "${CWB_SCRIPTS}/cwb_publish_ced.csh" 	#_SKIP_csh2sh_
#alias cwb_publish_postprod   "${CWB_SCRIPTS}/cwb_publish_postprod.csh"	#_SKIP_csh2sh_

alias cwb_report             "${CWB_SCRIPTS}/cwb_report.csh" 		#_SKIP_csh2sh_
alias cwb_mkhtml_prod        "${CWB_SCRIPTS}/cwb_mkhtml_prod.csh" 	#_SKIP_csh2sh_
alias cwb_mkhtml_sim         "${CWB_SCRIPTS}/cwb_mkhtml_sim.csh" 	#_SKIP_csh2sh_
alias cwb_mkeff              "${CWB_SCRIPTS}/cwb_mkeff.csh" 		#_SKIP_csh2sh_
alias cwb_cmpeff             "${CWB_SCRIPTS}/cwb_cmpeff.csh" 		#_SKIP_csh2sh_

alias cwb_report_cbc	     "${CWB_SCRIPTS}/cwb_report_cbc.csh"        #_SKIP_csh2sh_

alias cwb_setveto            "${CWB_SCRIPTS}/cwb_setveto.csh" 		#_SKIP_csh2sh_
alias cwb_setcuts            "${CWB_SCRIPTS}/cwb_setcuts.csh" 		#_SKIP_csh2sh_
alias cwb_combine_cbc        "${CWB_SCRIPTS}/cwb_combine_cbc.csh" 	#_SKIP_csh2sh_
alias cwb_setifar            "${CWB_SCRIPTS}/cwb_setifar.csh" 		#_SKIP_csh2sh_
alias cwb_setchunk           "${CWB_SCRIPTS}/cwb_setchunk.csh" 		#_SKIP_csh2sh_
alias cwb_setmulti           "${CWB_SCRIPTS}/cwb_setmulti.csh" 		#_SKIP_csh2sh_

alias cwb_csh2sh             "root -l -n ${CWB_MACROS}/cwb_csh2sh.C" 	#_SKIP_csh2sh_

alias cwb_dump               "${CWB_SCRIPTS}/cwb_dump.csh" 		#_SKIP_csh2sh_
alias cwb_draw_antpat        "${CWB_SCRIPTS}/cwb_draw_antpat.csh" 	#_SKIP_csh2sh_

alias cwb_draw_sensitivity   "${CWB_SCRIPTS}/cwb_draw_sensitivity.csh" 	#_SKIP_csh2sh_

alias cwb_fix_slag_missed    "${CWB_SCRIPTS}/cwb_fix_slag_missed.csh" 		#_SKIP_csh2sh_
alias cwb_fix_live_slag_missed "${CWB_SCRIPTS}/cwb_fix_live_slag_missed.csh" 	#_SKIP_csh2sh_

alias cwb_ls                 'ls ${CWB_MACROS}'				#_SKIP_csh2sh_
alias cwb_subcondor          'condor_submit_dag'			#_SKIP_csh2sh_

alias frdisplay              "${CWB_SCRIPTS}/cwb_frdisplay.csh" 	#_SKIP_csh2sh_

alias cwb_wavescan           "${CWB_WAVESCAN}/cwb_wavescan.csh" 	#_SKIP_csh2sh_

alias wavebico               "${HOME_WAVEBICO}/wavebico.csh" 		#_SKIP_csh2sh_

alias kroot  	             'ps | grep root | awk '"'"'{print $1}'"'"' | xargs kill -9'

alias cwb_mkfad              "${CWB_SCRIPTS}/cwb_mkfad.csh" 		#_SKIP_csh2sh_
alias cwb_mkrep              "${CWB_SCRIPTS}/cwb_mkrep.csh" 		#_SKIP_csh2sh_
alias cwb_mkhtml             "${CWB_SCRIPTS}/cwb_mkhtml.csh" 		#_SKIP_csh2sh_

alias cwb_help               "${CWB_SCRIPTS}/cwb_help.csh" 		#_SKIP_csh2sh_
alias cwb_watenv             "${CWB_SCRIPTS}/cwb_watenv.sh verbose" 

alias grid                   "grid-proxy-init"

if (! $?_USE_PEGASUS ) then
  alias cwb_condor           "${CWB_SCRIPTS}/cwb_condor.csh" 		#_SKIP_csh2sh_
  alias cwb_cstatus          "${CWB_SCRIPTS}/cwb_condor.csh status"	#_SKIP_csh2sh_
  alias cwb_ccreate          "${CWB_SCRIPTS}/cwb_condor.csh create"	#_SKIP_csh2sh_
  alias cwb_crecovery        "${CWB_SCRIPTS}/cwb_condor.csh recovery"	#_SKIP_csh2sh_
  alias cwb_csubmit          "${CWB_SCRIPTS}/cwb_condor.csh submit"	#_SKIP_csh2sh_
  alias cwb_cbench           "${CWB_SCRIPTS}/cwb_condor.csh benchmark"	#_SKIP_csh2sh_
  alias cwb_ccheck           "${CWB_SCRIPTS}/cwb_condor.csh check"	#_SKIP_csh2sh_
  alias cwb_clist            "${CWB_SCRIPTS}/cwb_condor.csh list"	#_SKIP_csh2sh_
  alias cwb_cremove          "${CWB_SCRIPTS}/cwb_condor.csh remove"	#_SKIP_csh2sh_
  alias cwb_ccleanup         "${CWB_SCRIPTS}/cwb_condor.csh cleanup"	#_SKIP_csh2sh_
else
  alias cwb_pegasus          "${CWB_SCRIPTS}/cwb_pegasus.csh"		#_SKIP_csh2sh_
  alias cwb_pcreate          "${CWB_SCRIPTS}/cwb_pegasus.csh create"	#_SKIP_csh2sh_
  alias cwb_precovery        "${CWB_SCRIPTS}/cwb_pegasus.csh recovery"	#_SKIP_csh2sh_
  alias cwb_psubmit          "${CWB_SCRIPTS}/cwb_pegasus.csh submit"	#_SKIP_csh2sh_
  alias cwb_pbench           "${CWB_SCRIPTS}/cwb_pegasus.csh benchmark"	#_SKIP_csh2sh_
  alias cwb_pstatus          "${CWB_SCRIPTS}/cwb_pegasus.csh status"	#_SKIP_csh2sh_
  alias cwb_premove          "${CWB_SCRIPTS}/cwb_pegasus.csh remove"	#_SKIP_csh2sh_
  alias cwb_prun             "${CWB_SCRIPTS}/cwb_pegasus.csh run"	#_SKIP_csh2sh_
  alias cwb_panalyzer        "${CWB_SCRIPTS}/cwb_pegasus.csh analyzer"	#_SKIP_csh2sh_
  alias cwb_pplots           "${CWB_SCRIPTS}/cwb_pegasus.csh plots"	#_SKIP_csh2sh_
  alias cwb_ptimeleft        "${CWB_SCRIPTS}/cwb_pegasus.csh timeleft" 	#_SKIP_csh2sh_
  alias voms                 "voms-proxy-init -voms virgo"
endif

if (! $?_USE_LSF ) then
  echo -n ""
else
  alias cwb_lsf              "${CWB_SCRIPTS}/cwb_lsf.csh"		#_SKIP_csh2sh_
endif

alias aladin 	             "$HOME_ALADIN/Aladin -nobanner"

source ${HOME_FRDISPLAY}/config.csh

unsetenv CWB_PARMS_FILES
unsetenv CWB_PPARMS_FILES
setenv CWB_PARMS_FILES	     '${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE}'
if (! $?_USE_ROOT6 ) then
# ROOT5 allows reinitialization of parameters
# we use both CWB_PARAMETERS_FILE,CWB_UPARAMETERS_FILE and CWB_MPARAMETERS_FILE for back-compatibility
# because CWB_MPARAMETERS_FILE has been introduced starting from SVN version
setenv CWB_PPARMS_FILES	     '${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_MPARAMETERS_FILE} ${CWB_EMPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE} ${CWB_PPARAMETERS_FILE} ${CWB_UPPARAMETERS_FILE} ${CWB_EPPARAMETERS_FILE}'
else
# CWB_MPARAMETERS_FILE contains reinitialization of parameters. 
# ROOT6 do not allows reinitialization of parameters, therefore
# CWB_MPARAMETERS_FILE can not be used
setenv CWB_PPARMS_FILES	     '${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_EMPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE} ${CWB_PPARAMETERS_FILE} ${CWB_UPPARAMETERS_FILE} ${CWB_EPPARAMETERS_FILE}'
endif

if (! $?ZSH_NAME ) then
  echo -n ""
else
  # this code is used by the zsh shell
  if ("$ZSH_NAME" == zsh) then
    eval "CWB_PARMS_FILES_ZSH=($CWB_PARMS_FILES)"
    eval "CWB_PPARMS_FILES_ZSH=($CWB_PPARMS_FILES)"
  endif
endif

