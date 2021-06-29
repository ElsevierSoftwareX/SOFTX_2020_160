#!/bin/bash

if [[ -z $SITE_CLUSTER  ]]; then
  echo "cwb setup error ..."
  echo "try 'source tools/site_config.sh' (where SITE_CLUSTER=ATLAS/CIT/CNAF/CASCINA/GATECH/USER_DEFINED)"
  return 0
fi

export ICC_INIT_SCRIPT="/opt/intel/2018/intel.sh"

if [[ "$SITE_CLUSTER" = ATLAS ]]; then

  if [[ -z $PYTHONPATH  ]]; then
    export PYTHONPATH=""
  fi
  if [[ -z $DYLD_LIBRARY_PATH  ]]; then
    export DYLD_LIBRARY_PATH=""
  fi
  if [[ "$DYLD_LIBRARY_PATH" != ${PYTHONPATH} ]]; then
    export DYLD_LIBRARY_PATH=${PYTHONPATH}:${DYLD_LIBRARY_PATH}
  fi
  if [[ ! $LD_LIBRARY_PATH =~ "${PYTHONPATH}" ]]; then
    export LD_LIBRARY_PATH=${PYTHONPATH}:${LD_LIBRARY_PATH}
  fi
  
  if [[ ! $PYTHONPATH =~ "${ROOTSYS}/lib:" ]]; then
    export PYTHONPATH=${ROOTSYS}/lib:${PYTHONPATH}
  fi
  if [[ ! $PYTHONPATH =~ "${HOME_LAL}/lib/python2.7/site-packages:" ]]; then
    export PYTHONPATH=${HOME_LAL}/lib/python2.7/site-packages:${PYTHONPATH}
  fi 
  if [[ ! $LD_LIBRARY_PATH =~ "/usr/lib/" ]]; then
    export LD_LIBRARY_PATH=/usr/lib/:${LD_LIBRARY_PATH}
  fi

  export NODE_DATA_DIR="/atlas/user/X_HOST_NAME/X_USER_NAME"
# export=NODE_DATA_DIR "/local/user/X_USER_NAME" # user node data name
# export=CONDOR_LOG_DIR "/local/user/X_USER_NAME" # condor log
  export CONDOR_LOG_DIR="/atlas/user/X_HOST_NAME/X_USER_NAME/condor"
  export WWW_PUBLIC_DIR="X_HOME/WWW/LSC/reports"

  export ICC_INIT_SCRIPT="/opt/intel/2018/intel.sh"

fi

if [[ "$SITE_CLUSTER" = CIT ]]; then

  if [[ -z $PYTHONPATH  ]]; then
    export PYTHONPATH=""
  fi
  if [[ -z $DYLD_LIBRARY_PATH  ]]; then
    export DYLD_LIBRARY_PATH=""
  fi
  if [[ "$DYLD_LIBRARY_PATH" != ${PYTHONPATH} ]]; then
    export DYLD_LIBRARY_PATH=${PYTHONPATH}:${DYLD_LIBRARY_PATH}
  fi
  if [[ ! $LD_LIBRARY_PATH =~ "${PYTHONPATH}" ]]; then
    export LD_LIBRARY_PATH=${PYTHONPATH}:${LD_LIBRARY_PATH}
  fi
  
  if [[ ! $PYTHONPATH =~ "${ROOTSYS}/lib:" ]]; then
    export PYTHONPATH=${ROOTSYS}/lib:${PYTHONPATH}
  fi
  if [[ ! $PYTHONPATH =~ "${HOME_LAL}/lib64/python2.7/site-packages:" ]]; then
    export PYTHONPATH=${HOME_LAL}/lib64/python2.7/site-packages:${PYTHONPATH}
  fi 
  if [[ ! $LD_LIBRARY_PATH =~ "/usr/lib/:/usr/lib64/" ]]; then
    export LD_LIBRARY_PATH=/usr/lib/:/usr/lib64/:${LD_LIBRARY_PATH}
  fi

# export=NODE_DATA_DIR "/data/X_HOST_NAME/X_USER_NAME" # user node data name
  export NODE_DATA_DIR="/local/X_USER_NAME"
  export CONDOR_LOG_DIR="/usr1/X_USER_NAME"
  export WWW_PUBLIC_DIR="X_HOME/public_html/reports"

  export ICC_INIT_SCRIPT="/ldcg/intel/2020u2/compilers_and_libraries_2020.2.254/linux/bin/iccvars.sh"

fi

if [[ "$SITE_CLUSTER" = CNAF ]]; then

  if [[ ! $LD_LIBRARY_PATH =~ "/usr/lib/:/usr/lib64/" ]]; then
    export LD_LIBRARY_PATH=/usr/lib/:/usr/lib64/:${LD_LIBRARY_PATH}
  fi

  export NODE_DATA_DIR="/home/VIRGO/X_USER_NAME"
  export CONDOR_LOG_DIR=""
  export WWW_PUBLIC_DIR=""

fi

if [[ "$SITE_CLUSTER" = CASCINA ]]; then

  if [[ ! $LD_LIBRARY_PATH =~ "/usr/lib/:/usr/lib64/" ]]; then
    export LD_LIBRARY_PATH=/usr/lib/:/usr/lib64/:${LD_LIBRARY_PATH}
  fi

  export NODE_DATA_DIR=""
  export CONDOR_LOG_DIR=""
  export WWW_PUBLIC_DIR=""

fi

if [[ "$SITE_CLUSTER" = GATECH ]]; then

  if [[ ! $LD_LIBRARY_PATH =~ "${METAIOLIB}" ]]; then
    export LD_LIBRARY_PATH=${METAIOLIB}:${LD_LIBRARY_PATH}
  fi

  export NODE_DATA_DIR="X_HOME"
  export CONDOR_LOG_DIR="condor"
  export WWW_PUBLIC_DIR=""

fi

if [[ "$SITE_CLUSTER" != CIT && "$SITE_CLUSTER" != ATLAS && "$SITE_CLUSTER" != CASCINA && "$SITE_CLUSTER" != CNAF && "$SITE_CLUSTER" != GATECH ]]; then

  export NODE_DATA_DIR=""
  export CONDOR_LOG_DIR=""
  export WWW_PUBLIC_DIR=""

fi

# define the batch system
if [[ -z $BATCH_SYSTEM  ]]; then
  export BATCH_SYSTEM=""
fi
if [[ "$BATCH_SYSTEM" != PEGASUS && "$BATCH_SYSTEM" != LSF && "$BATCH_SYSTEM" != OSG && "$BATCH_SYSTEM" != "" ]]; then
  export BATCH_SYSTEM=""
fi
if [[ "$BATCH_SYSTEM" = CONDOR || "$BATCH_SYSTEM" = "" ]]; then
  unset _USE_PEGASUS
  unset _USE_LSF
  unset _USE_OSG
fi
if [[ "$BATCH_SYSTEM" = PEGASUS ]]; then
  if [[ -z $CWB_PEGASUS_WATENV  ]]; then
    echo ""
    echo "env CWB_PEGASUS_WATENV is not defined !"
    return 0
  fi
  if [[ -z $CWB_PEGASUS_SITE  ]]; then
    echo ""
    echo "env CWB_PEGASUS_SITE is not defined !"
    return 0
  fi
  export _USE_PEGASUS=1
  unset _USE_LSF
  unset _USE_OSG
fi
if [[ "$BATCH_SYSTEM" = LSF ]]; then
  if [[ -z $LSF_QUEUE  ]]; then
    echo ""
    echo "env LSF_QUEUE is not defined !"
    return 0
  fi
  unset _USE_PEGASUS
  export _USE_LSF=1
  unset _USE_OSG
fi
if [[ "$BATCH_SYSTEM" = OSG ]]; then
  unset _USE_PEGASUS
  unset _USE_LSF
  export _USE_OSG=1
fi

if [[ -z $_USE_ICC  ]]; then
   echo -n ""
else
  if [[  ! -e $ICC_INIT_SCRIPT  ]]; then
    echo ""
    echo -n $ICC_INIT_SCRIPT
    echo " not exist!"
    echo "intel compiler can not be used, disable ICC in watenv -> unset _USE_ICC"
    return 0
  else
    if [[ ! $LD_LIBRARY_PATH =~ "intel64" ]]; then
      if [[ "$SITE_CLUSTER" = ATLAS ]]; then
        source $ICC_INIT_SCRIPT
      else 
        source $ICC_INIT_SCRIPT -arc intel64
      fi 
    fi 
#    if [[  $? != 0 ]]; then
#      echo ""
#      echo "source" $ICC_INIT_SCRIPT "intel64 error !!!"
#      echo "intel compiler can not be used, disable ICC in watenv -> unset _USE_ICC"
#      echo ""
#      return 0
#    fi
  fi
fi

if [[ -z $CWB_ANALYSIS  ]]; then
  echo ""
  echo "CWB_ANALYSIS env not defined in watenv setup"
  echo "add in watenv setup the following statement"
  echo ""
  echo "export=CWB_ANALYSIS '1G'
  echo "or"
  echo "export=CWB_ANALYSIS '2G'
  echo ""
  echo "use 'cwb_setpipe [1G/1g/1,2G/2g/2]' command to switch analysis type [1G/2G]"
  return 0
fi
if [[   ${CWB_ANALYSIS} != "1G"   &&   ${CWB_ANALYSIS} != "2G"   ]]; then
  echo -n "CWB_ANALYSIS="
  echo ${CWB_ANALYSIS}
  echo "CWB_ANALYSIS env must be 1G or 2G"
  echo "check watenv setup"
  return 0
fi

export CWB_TOOLBOX="${HOME_WAT}/tools/toolbox"
export CWB_SKYPLOT="${HOME_WAT}/tools/skyplot"
export CWB_GWAT="${HOME_WAT}/tools/gwat"
export CWB_EBBH="${HOME_WAT}/tools/eBBH"
export CWB_HISTORY="${HOME_WAT}/tools/history"
export CWB_STFT="${HOME_WAT}/tools/stft"
export CWB_BICO="${HOME_WAT}/tools/bico"
export CWB_FILTER="${HOME_WAT}/tools/filter"
export CWB_FRAME="${HOME_WAT}/tools/frame"
export CWB_WAVESCAN="${HOME_WAT}/tools/wavescan"
export CWB_ONLINE="${HOME_WAT}/tools/online"

if [[ ! $PATH =~ "${ROOTSYS}/bin" ]]; then
  export PATH="${ROOTSYS}/bin:${PATH}"
fi
if [[  "$OSTYPE" =~ "linux"  ]]; then
  if [[ ! $PATH =~ "${HOME_FRLIB}/Linux" ]]; then
    export PATH="${HOME_FRLIB}/Linux:${PATH}"
  fi
fi
if [[  "$OSTYPE" =~ "darwin"  ]]; then
  if [[ ! $PATH =~ "${HOME_FRLIB}/Darwin" ]]; then
    export PATH="${HOME_FRLIB}/Darwin:${PATH}"
  fi
fi
if [[ ! $LD_LIBRARY_PATH =~ "${ROOTSYS}/lib:" ]]; then
  export LD_LIBRARY_PATH="${ROOTSYS}/lib:${LD_LIBRARY_PATH}"
fi
if [[  $?_USE_ROOT6  ]]; then
  if [[ ! $LD_LIBRARY_PATH =~ "${HOME_WAT_INSTALL}/lib:" ]]; then
    export LD_LIBRARY_PATH="${HOME_WAT_INSTALL}/lib:${LD_LIBRARY_PATH}"
  fi
fi
if [[  $?_USE_LAL  ]]; then
  if [[ ! $LD_LIBRARY_PATH =~ "${HOME_LAL}/lib:" ]]; then
    export LD_LIBRARY_PATH="${HOME_LAL}/lib:${LD_LIBRARY_PATH}"
  fi
fi

# export=ROOT_INCLUDE_PATH: used for include path by macros
export ROOT_INCLUDE_PATH="$CWB_CONFIG"

if [[ -z $MANPATH  ]]; then
  export MANPATH=""
fi
if [[ ! $MANPATH =~ "${HOME_WAT}/man" ]]; then
  export MANPATH="${HOME_WAT}/man:${MANPATH}"
fi

alias root="${ROOTSYS}/bin/root -l"
alias broot="${ROOTSYS}/bin/root -b -l"

export HOME_CWB="${HOME_WAT}/tools/cwb"
export CWB_SCRIPTS="${HOME_CWB}/scripts"
export CWB_DIR_NAME="${HOME_CWB}/macros"
export CWB_MACROS="${HOME_CWB}/macros"

export HOME_WAVEBICO="${HOME_WAT}/tools/wavebico"

export HOME_WAVEMDC="${HOME_WAT}/tools/wavemdc"
export WMDC_CONFIG="./wmdc_config.C"
alias wmdc_condor="${HOME_WAVEMDC}/wmdc_condor.csh"

export HOME_FRDISPLAY="${HOME_WAT}/tools/frdisplay"

export CWB_ROOTLOGON_FILE="${CWB_MACROS}/cwb_rootlogon.C"

if [[  "$CWB_ANALYSIS" = "1G"  ]]; then
  export CWB_PARAMETERS_FILE="${CWB_MACROS}/cwb1G_parameters.C"
fi
if [[  "$CWB_ANALYSIS" = "2G"  ]]; then
  export CWB_PARAMETERS_FILE="${CWB_MACROS}/cwb2G_parameters.C"
fi
if [[ -z $CWB_UPARAMETERS_FILE  ]]; then
  export CWB_UPARAMETERS_FILE="./config/user_parameters.C"
fi
if [[ -z $CWB_MPARAMETERS_FILE  ]]; then
  export CWB_MPARAMETERS_FILE="./merge/cwb_user_parameters.C"
fi
if [[ -z $CWB_EMPARAMETERS_FILE  ]]; then
  export CWB_EMPARAMETERS_FILE=""
fi
export CWB_EPARAMETERS_FILE="${CWB_MACROS}/cwb_eparameters.C"
export CWB_PPARAMETERS_FILE="${CWB_MACROS}/cwb_pparameters.C"
if [[ -z $CWB_UPPARAMETERS_FILE  ]]; then
  export CWB_UPPARAMETERS_FILE="./config/user_pparameters.C"
fi
export CWB_EPPARAMETERS_FILE="${CWB_MACROS}/cwb_epparameters.C"

export CWB_NETS_FILE="${CWB_SCRIPTS}/cwb_net.sh"

if [[ -z $_USE_OSG  ]]; then
  echo -n ""
else
  export CWB_NETS_FILE="${CWB_SCRIPTS}/cwb_net_osg.sh"
fi
if [[ -z $_USE_LSF  ]]; then
  echo -n ""
else
  export CWB_NETS_FILE="${CWB_SCRIPTS}/cwb_net_lsf.sh"
fi
if [[ -z $_USE_PEGASUS  ]]; then
  echo -n ""
else
  export CWB_NETS_FILE="${CWB_SCRIPTS}/cwb_net_pegasus.sh"
fi

export CWB_NETC_FILE="${CWB_MACROS}/cwb_net.C"

export CWB_HTML_INDEX="${CWB_MACROS}/html_templates/html_index_template.txt"
export CWB_HTML_HEADER="${CWB_MACROS}/html_templates/html_header_template.txt"
export CWB_HTML_BODY_PROD="${CWB_MACROS}/html_templates/html_body_prod_template.txt"
export CWB_HTML_BODY_SIM="${CWB_MACROS}/html_templates/html_body_sim_template.txt"

export CWB_LAG_NUMBER=-1
export CWB_SLAG_NUMBER=-1

alias cwb_setpipe="source ${CWB_SCRIPTS}/cwb_setpipe.sh"
alias cwb_setpars="source ${CWB_SCRIPTS}/cwb_setpars.sh"
alias cwb_setppars="source ${CWB_SCRIPTS}/cwb_setppars.sh"
alias cwb_clonedir="${CWB_SCRIPTS}/cwb_clonedir.csh"

alias cwb_mkdir="${CWB_SCRIPTS}/cwb_mkdir.csh"
alias cwb_online="${CWB_SCRIPTS}/cwb_online.csh"
alias cwb_json2dat="${CWB_SCRIPTS}/cwb_json2dat.csh"
alias cwb_petools="${CWB_SCRIPTS}/cwb_petools.csh"
alias cwb_pereport="${CWB_SCRIPTS}/cwb_pereport.csh"

alias cwb_tune_slag="${CWB_SCRIPTS}/cwb_tune_slag.csh"
alias cwb_get_rejslag="${CWB_SCRIPTS}/cwb_get_rejslag.csh"

alias cwb_sub_lsf="${CWB_SCRIPTS}/cwb_sub_lsf.csh"

alias cwb_inet="${CWB_SCRIPTS}/cwb_inet.csh"
alias cwb_xnet="${CWB_SCRIPTS}/cwb_xnet.csh"
alias cwb_jnet="${CWB_SCRIPTS}/cwb_jnet.csh"
alias cwb_eced="${CWB_SCRIPTS}/cwb_eced.csh"
alias cwb_eced2G="${CWB_SCRIPTS}/cwb_eced2G.csh"
alias cwb_inet2G="${CWB_SCRIPTS}/cwb_inet2G.csh"
alias cwb_inet2G1="${CWB_SCRIPTS}/cwb_inet2G1.csh"
alias cwb_inet2G2="${CWB_SCRIPTS}/cwb_inet2G2.csh"
alias cwb_mplugin="${CWB_SCRIPTS}/cwb_mplugin.csh"
alias cwb_gwosc="make -f ${HOME_CWB}/gwosc/Makefile.gwosc"
alias cwb_compile="${CWB_SCRIPTS}/cwb_compile.csh"

alias cwb_merge="${CWB_SCRIPTS}/cwb_merge.csh"
alias cwb_xtalk="${CWB_SCRIPTS}/cwb_xtalk.csh"

alias cwb_publish="${CWB_SCRIPTS}/cwb_publish.csh"
alias cwb_publish_ced="${CWB_SCRIPTS}/cwb_publish_ced.csh"
#alias cwb_publish_postprod="${CWB_SCRIPTS}/cwb_publish_postprod.csh" #_SKIP_csh2sh_

alias cwb_report="${CWB_SCRIPTS}/cwb_report.csh"
alias cwb_mkhtml_prod="${CWB_SCRIPTS}/cwb_mkhtml_prod.csh"
alias cwb_mkhtml_sim="${CWB_SCRIPTS}/cwb_mkhtml_sim.csh"
alias cwb_mkeff="${CWB_SCRIPTS}/cwb_mkeff.csh"
alias cwb_cmpeff="${CWB_SCRIPTS}/cwb_cmpeff.csh"

alias cwb_report_cbc="${CWB_SCRIPTS}/cwb_report_cbc.csh"

alias cwb_setveto="${CWB_SCRIPTS}/cwb_setveto.csh"
alias cwb_setcuts="${CWB_SCRIPTS}/cwb_setcuts.csh"
alias cwb_combine_cbc="${CWB_SCRIPTS}/cwb_combine_cbc.csh"
alias cwb_setifar="${CWB_SCRIPTS}/cwb_setifar.csh"
alias cwb_setchunk="${CWB_SCRIPTS}/cwb_setchunk.csh"
alias cwb_setmulti="${CWB_SCRIPTS}/cwb_setmulti.csh"

alias cwb_csh2sh="root -l -n ${CWB_MACROS}/cwb_csh2sh.C"

alias cwb_dump="${CWB_SCRIPTS}/cwb_dump.csh"
alias cwb_draw_antpat="${CWB_SCRIPTS}/cwb_draw_antpat.csh"

alias cwb_draw_sensitivity="${CWB_SCRIPTS}/cwb_draw_sensitivity.csh"

alias cwb_fix_slag_missed="${CWB_SCRIPTS}/cwb_fix_slag_missed.csh"
alias cwb_fix_live_slag_missed="${CWB_SCRIPTS}/cwb_fix_live_slag_missed.csh"

alias cwb_ls='ls ${CWB_MACROS}'
alias cwb_subcondor='condor_submit_dag'

alias frdisplay="${CWB_SCRIPTS}/cwb_frdisplay.csh"

alias cwb_wavescan="${CWB_WAVESCAN}/cwb_wavescan.csh"

alias wavebico="${HOME_WAVEBICO}/wavebico.csh"

alias kroot='ps | grep root | awk '"'"'{print $1}'"'"' | xargs kill -9'

alias cwb_mkfad="${CWB_SCRIPTS}/cwb_mkfad.csh"
alias cwb_mkrep="${CWB_SCRIPTS}/cwb_mkrep.csh"
alias cwb_mkhtml="${CWB_SCRIPTS}/cwb_mkhtml.csh"

alias cwb_help="${CWB_SCRIPTS}/cwb_help.csh"
alias cwb_watenv="${CWB_SCRIPTS}/cwb_watenv.sh verbose"

alias grid="grid-proxy-init"

if [[ -z $_USE_PEGASUS  ]]; then
  alias cwb_condor="${CWB_SCRIPTS}/cwb_condor.csh"
  alias cwb_cstatus="${CWB_SCRIPTS}/cwb_condor.csh status"
  alias cwb_ccreate="${CWB_SCRIPTS}/cwb_condor.csh create"
  alias cwb_crecovery="${CWB_SCRIPTS}/cwb_condor.csh recovery"
  alias cwb_csubmit="${CWB_SCRIPTS}/cwb_condor.csh submit"
  alias cwb_cbench="${CWB_SCRIPTS}/cwb_condor.csh benchmark"
  alias cwb_ccheck="${CWB_SCRIPTS}/cwb_condor.csh check"
  alias cwb_clist="${CWB_SCRIPTS}/cwb_condor.csh list"
  alias cwb_cremove="${CWB_SCRIPTS}/cwb_condor.csh remove"
  alias cwb_ccleanup="${CWB_SCRIPTS}/cwb_condor.csh cleanup"
else
  alias cwb_pegasus="${CWB_SCRIPTS}/cwb_pegasus.csh"
  alias cwb_pcreate="${CWB_SCRIPTS}/cwb_pegasus.csh create"
  alias cwb_precovery="${CWB_SCRIPTS}/cwb_pegasus.csh recovery"
  alias cwb_psubmit="${CWB_SCRIPTS}/cwb_pegasus.csh submit"
  alias cwb_pbench="${CWB_SCRIPTS}/cwb_pegasus.csh benchmark"
  alias cwb_pstatus="${CWB_SCRIPTS}/cwb_pegasus.csh status"
  alias cwb_premove="${CWB_SCRIPTS}/cwb_pegasus.csh remove"
  alias cwb_prun="${CWB_SCRIPTS}/cwb_pegasus.csh run"
  alias cwb_panalyzer="${CWB_SCRIPTS}/cwb_pegasus.csh analyzer"
  alias cwb_pplots="${CWB_SCRIPTS}/cwb_pegasus.csh plots"
  alias cwb_ptimeleft="${CWB_SCRIPTS}/cwb_pegasus.csh timeleft"
  alias voms="voms-proxy-init -voms virgo"
fi

if [[ -z $_USE_LSF  ]]; then
  echo -n ""
else
  alias cwb_lsf="${CWB_SCRIPTS}/cwb_lsf.csh"
fi

alias aladin="$HOME_ALADIN/Aladin -nobanner"

source ${HOME_FRDISPLAY}/config.sh

unset CWB_PARMS_FILES
unset CWB_PPARMS_FILES
export CWB_PARMS_FILES='${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE}'
if [[ -z $_USE_ROOT6  ]]; then
# ROOT5 allows reinitialization of parameters
# we use both CWB_PARAMETERS_FILE,CWB_UPARAMETERS_FILE and CWB_MPARAMETERS_FILE for back-compatibility
# because CWB_MPARAMETERS_FILE has been introduced starting from SVN version
export CWB_PPARMS_FILES='${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_MPARAMETERS_FILE} ${CWB_EMPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE} ${CWB_PPARAMETERS_FILE} ${CWB_UPPARAMETERS_FILE} ${CWB_EPPARAMETERS_FILE}'
else
# CWB_MPARAMETERS_FILE contains reinitialization of parameters. 
# ROOT6 do not allows reinitialization of parameters, therefore
# CWB_MPARAMETERS_FILE can not be used
export CWB_PPARMS_FILES='${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_EMPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE} ${CWB_PPARAMETERS_FILE} ${CWB_UPPARAMETERS_FILE} ${CWB_EPPARAMETERS_FILE}'
fi

if [[ -z $ZSH_NAME  ]]; then
  echo -n ""
else
  
  if [[ "$ZSH_NAME" = zsh ]]; then
    eval "CWB_PARMS_FILES_ZSH=($CWB_PARMS_FILES)"
    eval "CWB_PPARMS_FILES_ZSH=($CWB_PPARMS_FILES)"
  fi
fi

