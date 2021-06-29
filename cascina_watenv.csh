#  -----------------------------------------------------------------
#  WAT env for CASCINA CLUSTER
#  -----------------------------------------------------------------

  if (! $?LD_LIBRARY_PATH ) then
    setenv LD_LIBRARY_PATH ""
  endif

  setenv SITE_CLUSTER	     "CASCINA"
  setenv BATCH_SYSTEM        "CONDOR"

  unsetenv _USE_ICC          # disable ICC
  #setenv   _USE_ICC 1        # enable ICC   
  #unsetenv _USE_CPP11        # disable C++11
  setenv _USE_CPP11 1        # enable C++11 
  #unsetenv _USE_ROOT6        # disable ROOT6
  setenv _USE_ROOT6 1        # enable ROOT6 (gcc must be > 4.8.2, enable C++11)

#  -----------------------------------------------------------------
#  this section is mandatory to compile the WAT libraries
#  -----------------------------------------------------------------

  # LIBRARIES PATH
  setenv HOME_LIBS           "/virgoDev/waveburst"

  # ROOT
  setenv ROOT_VERSION        "root-v6-14-06"
  setenv ROOTSYS 	     "/virgoApp/root/v6r14p060/Linux-x86_64-CL7"

  # FRAMELIB
  setenv HOME_FRLIB          "/virgo/ext/libframe/libframe-8.33_root-6.14.06_icc"

  # HEALPIX
  #unsetenv _USE_HEALPIX      # disable HEALPix in CWB
  setenv _USE_HEALPIX 1      # enable HEALPix in CWB
  setenv HOME_HEALPIX        "/virgo/ext/Healpix/Healpix_3.40/Linux-x86_64-CL7/"
  setenv HOME_CFITSIO        "/virgo/ext/cfitsio/cfitsio-3.45/Linux-x86_64-CL7/lib"

  # LAL
  unsetenv _USE_LAL          # disable LAL in CWB
  #setenv _USE_LAL 1          # enable LAL in CWB
  setenv HOME_LAL            "${HOME_LIBS}/LAL/lalsuite/master"
  setenv LAL_INC             "${HOME_LAL}/include"
  setenv LALINSPINJ_EXEC     "${HOME_LAL}/bin/lalapps_inspinj"

  # EBBH/CVODE
  unsetenv _USE_EBBH          # disable EBBH in CWB
  #setenv _USE_EBBH 1          # enable EBBH in CWB
  setenv HOME_CVODE            "${HOME_LIBS}/CVODE/cvode-2.7.0/dist"

#  -----------------------------------------------------------------
#  this section is specific for the CWB pipeline 
#  -----------------------------------------------------------------

  #setenv CWB_ANALYSIS        "1G"       # 1G  analysis
  setenv CWB_ANALYSIS        "2G"      # 2G  analysis
  
  setenv CWB_CONFIG          "CWB_CONFIG ENV NOT DEFINED!!! "   # Location of CWB Configuration files, DQ files, FRAMES files

  setenv HOME_WWW            "~waveburst/waveburst/WWW"
  setenv HOME_CED_WWW        "~waveburst/waveburst/WWW/ced"
  setenv HOME_CED_PATH       "/opt/w3/DataAnalysis/Burst/cWB/ced"
  setenv CWB_DOC_URL         "https://ldas-jobs.ligo.caltech.edu/~waveburst/waveburst/WWW/doc"
  setenv CWB_REP_URL         "https://ldas-jobs.ligo.caltech.edu/~waveburst/waveburst/WWW/reports"
  setenv CWB_GIT_URL         "https://git.ligo.org/cWB"

  setenv HOME_WAT_FILTERS    "/virgoApp/cWB/v6r2p6/WAT/filters"
  setenv HOME_BAUDLINE       "${HOME_LIBS}/BAUDLINE/baudline_1.08_linux_x86_64"
  setenv HOME_ALADIN         "${HOME_LIBS}/ALADIN/Aladin.v7.533"

  # select watenv script to be executed remote site : used when _USE_PEGASUS is enabled
  setenv CWB_PEGASUS_WATENV  "/cvmfs/virgo.infn.it/waveburst/SOFT/WAT/trunk/cnaf_watenv.sh"
  setenv CWB_PEGASUS_SITE    "creamce_cnaf"

#  -----------------------------------------------------------------
#  DO NOT MODIFY !!!
#
#  1) In this section the HOME_WAT env is automatically initialized 
#  2) Init Default cWB library 
#  3) Init Default cWB config  
#  -----------------------------------------------------------------

  # HOME_WAT is initalized with directory path name of the this script
  set MYSHELL=`readlink /proc/$$/exe`
  if ( "$MYSHELL" =~ "*tcsh" ) then
    echo "\nEntering in TCSH section..."
    set CWB_CALLED=($_)
    if ( $#CWB_CALLED == 2 ) then       # alias
      set CWB_CALLED=`alias $CWB_CALLED`
    endif
    if ( "$CWB_CALLED" != "" ) then     # called by source 
      set WATENV_SCRIPT=`readlink -f $CWB_CALLED[2]`
    else                                # called by direct excution of the script
      echo "\nError: script must be executed with source command"
      exit 1
    endif
    set script_dir=`dirname $WATENV_SCRIPT`
  endif
  if ( "$MYSHELL" =~ "*bash" ) then
    echo ""
    echo "Entering in BASH section..."
    script_dir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ) #_SKIP_CSH2SH_
  endif
  if ( "$MYSHELL" =~ "*zsh" ) then
    echo "\nEntering in ZSH section..."
    script_dir=$( cd "$( dirname "${(%):-%N}" )" && pwd )        #_SKIP_CSH2SH_
  endif
  setenv HOME_WAT           $script_dir
  setenv HOME_WAT_INSTALL   "${HOME_WAT}/tools/install" # define the wat lib/inc
  echo ""
  echo "------------------------------------------------------------------------"
  echo " -> HOME_WAT = $HOME_WAT"
  echo "------------------------------------------------------------------------"
  echo ""

  source $HOME_WAT/tools/config.csh     # init default cWB library
#  source $CWB_CONFIG/setup.csh          # init default cWB config


