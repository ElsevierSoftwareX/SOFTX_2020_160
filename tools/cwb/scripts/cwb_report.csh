#!/bin/tcsh -f

onintr irq_ctrlc

set cmd_line="$0 $argv"

if ($1 == '') then
  $CWB_SCRIPTS/cwb_help.csh cwb_report  
  exit
endif

if ($2 == 'cbc') then
  ${CWB_SCRIPTS}/cwb_report_cbc.csh $1 $3 $4
  exit
endif

if ($1 == 'slags') then

  setenv CWB_MERGE_LABEL $2

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

    (${CWB_XPARM};${CWB_XPPARM};root -b -n -l ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_report_slags.C)
    if ( $? != 0 ) then
      echo "cwb_report_slags.C error : process terminated"
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

  exit
endif

if ($1 == 'skymap') then

  if ( "$2" != '' ) then
    setenv CWB_SKYMAP_FILE "$2"
  else
    unsetenv CWB_SKYMAP_FILE
    echo ""
    echo "cwb_report skymap error : no input skymap file, process terminated"
    echo ""
    exit 
  endif

  if ( "$3" != '' ) then
    setenv CWB_REPORT_OPTIONS "$3"
  else
    unsetenv CWB_REPORT_OPTIONS
  endif

  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_report_skymap.C
  if ( $? != 0) then
    echo ""
    echo "cwb_report_skymap.C error : process terminated"
    echo ""
    exit 
  endif

  # create cWB_analysis.log file
  make -f $CWB_SCRIPTS/Makefile.log CMD_LINE="$cmd_line" svn >& /dev/null

  unsetenv CWB_REPORT_OPTIONS
  unsetenv CWB_SKYMAP_FILE
  exit
endif

if ($1 == 'list') then
  if ($2 == 'merge') then
    root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_merge_dir.C
  else 
    root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_report_dir.C
  endif 
  exit
endif

setenv CWB_MERGE_LABEL $1
unsetenv CWB_REPORT_PE 

if ($2 == 'loudest') then

  if ( "$3" != '' ) then
    setenv CWB_REPORT_OPTIONS "$3"
  else
    unsetenv CWB_REPORT_OPTIONS
  endif

  setenv CWB_LAG_NUMBER -1
  setenv CWB_SLAG_NUMBER -1

  root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_report_loudest.C
  if ( $? != 0) then
    echo ""
    echo "cwb_report_loudest.C error : process terminated"
    echo ""
    exit 
  endif

  # create cWB_analysis.log file
  make -f $CWB_SCRIPTS/Makefile.log CMD_LINE="$cmd_line" svn >& /dev/null

  unsetenv CWB_REPORT_OPTIONS
  unsetenv CWB_MERGE_LABEL
  unsetenv CWB_SLAG_NUMBER
  unsetenv CWB_LAG_NUMBER
  exit
endif

if ($2 == 'remove') then
  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_rm_report_dir.C
  if ( $? != 0) then
    echo ""
    echo "cwb_rm_report_dir.C error : process terminated"
    echo ""
    exit 
  endif
  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_report_dir.C
  exit
endif

if (($2 == 'publish') || ($2 == 'clean')) then
  setenv CWB_PP_DATA_DIR $1
  setenv CWB_PUBLISH_OPTION $2
  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_publish.C
  if ( $? != 0) then
    echo ""
    echo "cwb_publish.C error : process terminated"
    echo ""
    exit 
  endif
  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_dump_report_dir.C
  exit
endif

if (($2 != 'create') && ($2 != 'pe')) then
  echo ""
  echo \'$1\' "is a wrong cwb_report option"
  echo "type cwb_report & return to see the availables options"
  echo ""
  exit
endif

if (($3 != '') && ( $3 =~ *[^0-9]* )) then
  setenv CWB_XPARAMETERS_FILE $3
else
  setenv CWB_XPARAMETERS_FILE ""
endif

if (($4 != '') && ( $4 =~ *[^0-9]* )) then
  setenv CWB_XPPARAMETERS_FILE $4
else
  setenv CWB_XPPARAMETERS_FILE ${CWB_UPPARAMETERS_FILE}
endif

set CWB_XPARM  = "setenv CWB_EMPARAMETERS_FILE ${CWB_XPARAMETERS_FILE}"
set CWB_XPPARM = "setenv CWB_UPPARAMETERS_FILE ${CWB_XPPARAMETERS_FILE}"


# pe

root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_getsim.C
if (( $? == 0 ) && ($2 == 'pe')) then
  echo ""
  echo "cwb_report error : pe option can be used only with simulation"
  echo ""
  exit
endif

root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_getsim.C
if (( $? == 1 ) && ($2 == 'pe')) then
  setenv CWB_LAG_NUMBER -1
  setenv CWB_SLAG_NUMBER -1
  setenv CWB_REPORT_PE 1

  (${CWB_XPARM};${CWB_XPPARM};root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_report_pe.C)
  if ( $? != 0) then
    echo ""
    echo "cwb_report_pe.C error : process terminated"
    echo ""
    exit 
  endif
  unsetenv CWB_MKREP_INDEX_FILE
  (${CWB_XPARM};${CWB_XPPARM};root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_mkhtml_index.C)
  if ( $? != 0) then
    echo ""
    echo "cwb_mkhtml_index.C error : process terminated"
    echo ""
    exit 
  endif
  (${CWB_XPARM};${CWB_XPPARM};root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_mkhtml_header.C)
  if ( $? != 0) then
    echo ""
    echo "cwb_mkhtml_header.C error : process terminated"
    echo ""
    exit 
  endif
  (${CWB_XPARM};${CWB_XPPARM};root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_mkhtml_pe.C)
  if ( $? != 0) then
    echo "cwb_mkhtml_pe.C error : process terminated"
    exit 
  endif

  # create cWB_analysis.log file
  make -f $CWB_SCRIPTS/Makefile.log CMD_LINE="$cmd_line" svn >& /dev/null

  unsetenv CWB_REPORT_PE 
  unsetenv CWB_MERGE_LABEL
  unsetenv CWB_SLAG_NUMBER
  unsetenv CWB_LAG_NUMBER
  unsetenv CWB_XPARAMETERS_FILE
  unsetenv CWB_XPPARAMETERS_FILE
  exit
endif

# simulation

root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_getsim.C
if ( $? == 1 ) then
  setenv CWB_LAG_NUMBER -1
  setenv CWB_SLAG_NUMBER -1

  (${CWB_XPARM};${CWB_XPPARM};root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_report_sim.C)
  if ( $? != 0) then
    echo ""
    echo "cwb_report_sim.C error : process terminated"
    echo ""
    exit 
  endif
  (${CWB_XPARM};${CWB_XPPARM};root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_mkeff.C)
  if ( $? != 0) then
    echo ""
    echo "cwb_mkeff.C error : process terminated"
    echo ""
    exit 
  endif
  unsetenv CWB_MKREP_INDEX_FILE
  (${CWB_XPARM};${CWB_XPPARM};root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_mkhtml_index.C)
  if ( $? != 0) then
    echo ""
    echo "cwb_mkhtml_index.C error : process terminated"
    echo ""
    exit 
  endif
  (${CWB_XPARM};${CWB_XPPARM};root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_mkhtml_header.C)
  if ( $? != 0) then
    echo ""
    echo "cwb_mkhtml_header.C error : process terminated"
    echo ""
    exit 
  endif
  (${CWB_XPARM};${CWB_XPPARM};root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_mkhtml_sim.C)
  if ( $? != 0) then
    echo "cwb_mkhtml_sim.C error : process terminated"
    exit 
  endif

  # create cWB_analysis.log file
  make -f $CWB_SCRIPTS/Makefile.log CMD_LINE="$cmd_line" svn >& /dev/null

  unsetenv CWB_MERGE_LABEL
  unsetenv CWB_SLAG_NUMBER
  unsetenv CWB_LAG_NUMBER
  unsetenv CWB_XPARAMETERS_FILE
  unsetenv CWB_XPPARAMETERS_FILE
  exit
endif

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

  (${CWB_XPARM};${CWB_XPPARM};root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_report_prod_1.C)
  if ( $? != 0) then
    echo ""
    echo "cwb_report_prod_1.C error : process terminated"
    echo ""
    exit 
  endif
  (${CWB_XPARM};${CWB_XPPARM};root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_report_prod_2.C)
  if ( $? != 0) then
    echo "cwb_report_prod_2.C error : process terminated"
    exit 
  endif
  unsetenv CWB_MKREP_INDEX_FILE
  (${CWB_XPARM};${CWB_XPPARM};root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_mkhtml_index.C)
  if ( $? != 0) then
    echo ""
    echo "cwb_mkhtml_index.C error : process terminated"
    echo ""
    exit 
  endif
  (${CWB_XPARM};${CWB_XPPARM};root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_mkhtml_header.C)
  if ( $? != 0) then
    echo ""
    echo "cwb_mkhtml_header.C error : process terminated"
    echo ""
    exit 
  endif
  (${CWB_XPARM};${CWB_XPPARM};root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_mkhtml_prod.C)
  if ( $? != 0) then
    echo ""
    echo "cwb_mkhtml_prod.C error : process terminated"
    echo ""
    exit 
  endif
  (${CWB_XPARM};${CWB_XPPARM};root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_condor_create_ced.C)
  if ( $? != 0) then
    echo ""
    echo "cwb_condor_create_ced.C error : process terminated"
    echo ""
    exit
  endif

  # create cWB_analysis.log file
  make -f $CWB_SCRIPTS/Makefile.log CMD_LINE="$cmd_line" svn >& /dev/null

  unsetenv CWB_MERGE_LABEL
  unsetenv CWB_SLAG_NUMBER
  unsetenv CWB_LAG_NUMBER
  unsetenv CWB_XPARAMETERS_FILE
  unsetenv CWB_XPPARAMETERS_FILE
  exit
endif

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

