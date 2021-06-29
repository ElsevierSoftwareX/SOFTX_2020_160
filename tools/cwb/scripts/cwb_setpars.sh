#!/bin/bash -f

if [[ -z $1 ]]; then
  $CWB_SCRIPTS/cwb_help.sh cwb_setpars
  return 0
fi

if [[ $1 = "q" ]]; then
  echo ""
  echo "Current CWB_UPARAMETERS_FILE :"
  echo -n "-> "
  echo ${CWB_UPARAMETERS_FILE}
  echo ""
  return 0
fi

if [[ $1 = "." ]]; then
  export CWB_UPARAMETERS_FILE=config/user_parameters.C
else
  export CWB_UPARAMETERS_FILE=$1
fi

echo ""
echo "Current CWB_UPARAMETERS_FILE :"
echo -n "-> "
echo ${CWB_UPARAMETERS_FILE}
echo ""
