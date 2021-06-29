#  -----------------------------------------------------------------
#  WAT env for CNAF CLUSTER
#  -----------------------------------------------------------------

  if [[ -z $LD_LIBRARY_PATH  ]]; then
    export LD_LIBRARY_PATH=""
  fi

  export SITE_CLUSTER="CNAF"
  export BATCH_SYSTEM="LSF"

  
  export _USE_ICC=1
  
  export _USE_CPP11=1
  
  export _USE_ROOT6=1
  unset _USE_PEGASUS
  

#  -----------------------------------------------------------------
#  this section is mandatory to compile the WAT libraries
#  -----------------------------------------------------------------

  
  export HOME_LIBS="/opt/exp_software/virgo/virgoDev/waveburst/SOFT"

  
  export ROOT_VERSION="root_v6.14.06"
  
  export ROOTSYS="${HOME_LIBS}/ROOT/${ROOT_VERSION}"

  
  export HOME_FRLIB="${HOME_LIBS}/FRAMELIB/libframe-8.33_root-6.14.06"
  

  
  export _USE_HEALPIX=1
  export HOME_HEALPIX="${HOME_LIBS}/HEALPix/Healpix_3.40"
  export HOME_CFITSIO="${HOME_LIBS}/CFITSIO/cfitsio-3.45"

  
  unset _USE_LAL          
  
  export HOME_LAL="${HOME_LIBS}/LAL/lalsuite_lal-v6.19.2"
  export LAL_INC="${HOME_LAL}/include"
  export LALINSPINJ_EXEC="${HOME_LAL}/bin/lalapps_inspinj"

  
  unset _USE_EBBH          
  
  export HOME_CVODE="${HOME_LIBS}/CVODE/sundials-2.7.0/dist"

#  -----------------------------------------------------------------
#  this section is specific for the CWB pipeline 
#  -----------------------------------------------------------------

  
  export CWB_ANALYSIS="2G"
  
  export CWB_CONFIG="${HOME_LIBS}/cWB/tags/config/cWB-config-cWB-6.3.0"

  export HOME_WWW="~waveburst/waveburst/WWW"
  export HOME_CED_WWW="~waveburst/waveburst/WWW/ced"
  export HOME_CED_PATH="/home/waveburst/public_html/waveburst/WWW/ced"
  export CWB_DOC_URL="https://ldas-jobs.ligo.caltech.edu/~waveburst/waveburst/WWW/doc"
  export CWB_REP_URL="https://ldas-jobs.ligo.caltech.edu/~waveburst/waveburst/WWW/reports"
  export CWB_GIT_URL="https://git.ligo.org/cWB"

  export HOME_WAT_FILTERS="${HOME_LIBS}/WAT/filters"
  export HOME_BAUDLINE="${HOME_LIBS}/BAUDLINE/baudline_1.08_linux_x86_64"
  export HOME_ALADIN="${HOME_LIBS}/ALADIN/Aladin.v7.533"

  
  export CWB_PEGASUS_WATENV="/cvmfs/virgo.infn.it/waveburst/SOFT/WAT/trunk/cnaf_watenv.sh"
  export CWB_PEGASUS_SITE="creamce_cnaf"

  
  export LSF_QUEUE="virgo"

#  -----------------------------------------------------------------
#  DO NOT MODIFY !!!
#
#  1) In this section the HOME_WAT env is automatically initialized 
#  2) Init Default cWB library 
#  3) Init Default cWB config  
#  -----------------------------------------------------------------

  
   MYSHELL=`readlink /proc/$$/exe`
  if [[  "$MYSHELL" =~ "tcsh"  ]]; then
    echo "\nEntering in TCSH section..."
     CWB_CALLED=($_)
    if [[=$#CWB_CALLED = 2 ]]; then # alias
      set CWB_CALLED=`alias=$CWB_CALLED`
    fi
    if [[  "$CWB_CALLED" != ""  ]]; then     
       WATENV_SCRIPT=`readlink -f $CWB_CALLED[2]`
    else                                
      echo "\nError: script must be executed with source command"
      return 0 1
    fi
     script_dir=`dirname $WATENV_SCRIPT`
  fi
  if [[  "$MYSHELL" =~ "bash"  ]]; then
    echo ""
    echo "Entering in BASH section..."
    script_dir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ) #_SKIP_CSH2SH_
  fi
  if [[  "$MYSHELL" =~ "zsh"  ]]; then
    echo "\nEntering in ZSH section..."
    script_dir=$( cd "$( dirname "${(%):-%N}" )" && pwd )        #_SKIP_CSH2SH_
  fi
  export HOME_WAT=$script_dir
  export HOME_WAT_INSTALL="${HOME_WAT}/tools/install"

  source $HOME_WAT/tools/config.sh     
  
  source $CWB_SCRIPTS/cwb_watenv.sh     

  
  export CWB_HTML_INDEX="${CWB_MACROS}/html_templates/html_index_template_modern.txt"

