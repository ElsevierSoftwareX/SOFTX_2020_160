#!/bin/tcsh -f

onintr irq_ctrlc

setenv CWB_MERGE_LABEL $1

if (($3 != '') && ( $3 =~ *[^0-9]* )) then
  setenv CWB_XPARAMETERS_FILE $3
else
  setenv CWB_XPARAMETERS_FILE ${CWB_UPARAMETERS_FILE}
endif

if (($4 != '') && ( $4 =~ *[^0-9]* )) then
  setenv CWB_XPPARAMETERS_FILE $4
else
  setenv CWB_XPPARAMETERS_FILE ${CWB_UPPARAMETERS_FILE}
endif

set CWB_XPARM  = "setenv CWB_UPARAMETERS_FILE  ${CWB_XPARAMETERS_FILE}"
set CWB_XPPARM = "setenv CWB_UPPARAMETERS_FILE ${CWB_XPPARAMETERS_FILE}"

root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_getsim.C


#production

if ( $? == 0 ) then
  if ($3 == '') then
    setenv CWB_LAG_NUMBER -1
  else  
    setenv CWB_LAG_NUMBER $3
  endif

  if ($4 == '') then
    setenv CWB_SLAG_NUMBER -1
  else  
    setenv CWB_SLAG_NUMBER $4
  endif

  #(${CWB_XPARM};${CWB_XPPARM};root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/../postproduction/burst/slag.C)
  (${CWB_XPARM};${CWB_XPPARM};root -n -l ${CWB_PPARMS_FILES} ${CWB_MACROS}/../postproduction/burst/slag.C)
  if ( $? != 0 ) then
    echo "slag.C error : process terminated"
    exit
  endif
  unsetenv CWB_MKREP_INDEX_FILE

  unsetenv CWB_MERGE_LABEL
  unsetenv CWB_SLAG_NUMBER
  unsetenv CWB_LAG_NUMBER
  unsetenv CWB_XPARAMETERS_FILE
  unsetenv CWB_XPPARAMETERS_FILE
  exit
endif

exit 0
irq_ctrlc:
  ps r | grep root | awk '{print $1}' | xargs kill -9
  exit 1

