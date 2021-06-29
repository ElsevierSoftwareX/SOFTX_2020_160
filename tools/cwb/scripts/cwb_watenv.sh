#!/bin/bash -f

# -------------------------------------------------------------------------------------------
# PRINT cWB ENVIRONMENT VARIABLES
# -------------------------------------------------------------------------------------------

  echo ""
  echo "************************************************************************************"
  echo "                              cWB Environment Variables                             "
  echo "************************************************************************************"
  echo ""
  echo " -> SITE             = $SITE_CLUSTER"
  echo "------------------------------------------------------------------------------------"
if [ -d $HOME_WAT/.git ]; then
  echo " -> HOME_WAT         = $HOME_WAT (`git -C $HOME_WAT branch | grep \* | cut -d ' ' -f2`)"
else
  echo " -> HOME_WAT         = $HOME_WAT"
fi
if [ -d $CWB_CONFIG/.git ]; then
  echo " -> CWB_CONFIG       = $CWB_CONFIG (`git -C $CWB_CONFIG branch | grep \* | cut -d ' ' -f2`)"
else
  echo " -> CWB_CONFIG       = $CWB_CONFIG"
fi
  echo "------------------------------------------------------------------------------------"
if [[ "$1" = "verbose" ]]; then
  echo " -> HOME_LIBS        = $HOME_LIBS"
  echo " -> ROOTSYS          = $ROOTSYS"
  echo " -> HOME_LAL         = $HOME_LAL"
  echo " -> HOME_FRLIB       = $HOME_FRLIB"
  echo " -> HOME_CFITSIO     = $HOME_CFITSIO"
  echo " -> HOME_HEALPIX     = $HOME_HEALPIX"
  echo " -> HOME_CVODE       = $HOME_CVODE"
  echo " -> HOME_WAT_FILTERS = $HOME_WAT_FILTERS"
  echo " -> HOME_BAUDLINE    = $HOME_BAUDLINE"
  echo " -> HOME_ALADIN      = $HOME_ALADIN"
  echo " -> HOME_SKYMAP_LIB  = $HOME_SKYMAP_LIB"
  echo "------------------------------------------------------------------------------------"
  echo " -> HOME_WWW         = $HOME_WWW"
  echo " -> CWB_USER_URL     = $CWB_USER_URL"
  echo " -> HOME_CED_WWW     = $HOME_CED_WWW"
  echo " -> HOME_CED_PATH    = $HOME_CED_PATH"
  echo " -> CWB_REP_URL      = $CWB_REP_URL"
  echo " -> CWB_DOC_URL      = $CWB_DOC_URL"
  echo " -> CWB_GIT_URL      = $CWB_GIT_URL"
  echo "------------------------------------------------------------------------------------"
if [[ $_USE_ROOT6  ]]; then
  echo " -> ROOT6              ENABLED"
else
  echo " -> ROOT6              DISABLE"
fi
if [[ $_USE_LAL  ]]; then
  echo " -> LAL                ENABLED"
else
  echo " -> LAL                DISABLE"
fi
if [[ $_USE_EBBH  ]]; then
  echo " -> EBBH               ENABLED"
else
  echo " -> EBBH               DISABLE"
fi
if [[ $_USE_HEALPIX  ]]; then
  echo " -> HEALPIX            ENABLED"
else
  echo " -> HEALPIX            DISABLE"
fi
if [[ $_USE_ICC  ]]; then
  echo " -> ICC                ENABLED"
else
  echo " -> ICC                DISABLE"
fi
  echo "------------------------------------------------------------------------------------"

fi
  echo ""

