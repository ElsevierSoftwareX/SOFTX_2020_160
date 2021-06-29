#!/bin/bash -f

if [[  $1 != "1G"  &&  $1 != "1g"  &&  $1 != "1"  &&  $1 != "2G"  &&  $1 != "2g"  &&  $1 != "2"  ]]; then
  $CWB_SCRIPTS/cwb_help.sh cwb_setpipe
  echo ""
  echo -n "Current CWB_ANALYSIS="
  echo ${CWB_ANALYSIS}
  echo ""
  return 0
fi

if [[  $1 = "1G"  ||  $1 = "1g"  ||  $1 = "1"  ]]; then
  export CWB_ANALYSIS="1G"
  source $HOME_WAT/tools/config.sh
  echo ""
  echo -n "CWB_ANALYSIS="
  echo ${CWB_ANALYSIS}
  echo ""
fi

if [[  $1 = "2G"  ||  $1 = "2g"  ||  $1 = "2"  ]]; then
  export CWB_ANALYSIS="2G"
  source $HOME_WAT/tools/config.sh
  echo ""
  echo -n "CWB_ANALYSIS="
  echo ${CWB_ANALYSIS}
  echo ""
fi

