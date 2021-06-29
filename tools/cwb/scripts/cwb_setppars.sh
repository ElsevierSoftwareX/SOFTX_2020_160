#!/bin/bash -f

if [[ -z $1 ]]; then
  $CWB_SCRIPTS/cwb_help.sh cwb_setppars
  return 0
fi

if [[ $1 = "q" ]]; then
  echo ""
  echo "Current CWB_UPPARAMETERS_FILE :"
  echo -n "-> "
  echo ${CWB_UPPARAMETERS_FILE}
  echo ""
  return 0
fi

if [[ $1 = "." ]]; then
  export CWB_UPPARAMETERS_FILE=config/user_pparameters.C
else
  export CWB_UPPARAMETERS_FILE=$1
fi

echo ""
echo "Current CWB_UPPARAMETERS_FILE :"
echo -n "-> "
echo ${CWB_UPPARAMETERS_FILE}
echo ""
