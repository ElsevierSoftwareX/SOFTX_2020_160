#!/bin/tcsh -f

onintr irq_ctrlc

set cmd_line="$0 $argv"

if ($1 == '') then
  $CWB_SCRIPTS/cwb_help.csh cwb_condor
  exit
endif

if ((( $1 != 'create' ) && ( $1 != 'submit' ) && ( $1 != 'recovery' ) && ( $1 != 'status' ) && ( $1 != 'check' ) && ( $1 != 'list' ) && ( $1 != 'benchmark' ) && ( $1 != 'cleanup' ) && ( $1 != 'remove' ) && ( $1 != 'resume' ) && ( $1 != 'sdag' )  && ( $1 != 'xdag' ) && ( $1 != 'mtpe' ))) then
  echo ""
  echo \'$1\' "is a wrong cwb_condor option"
  echo "type cwb_condor & return to see the availables options"
  echo ""
  exit
endif

if ($1 == 'mtpe') then  
  if ($2 == '') then
    echo ""
    echo "cwb_condor_mtpe - Error : missing parameter"
    echo "Example : cwb_condor mtpe JOBID"
    echo "          where JOBD is the JOB NUMBER"
    echo ""
    exit
  else
    setenv CWB_CONDOR_MTPE_JOBID $2
    root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_condor_mtpe.C
    unsetenv CWB_CONDOR_MTPE_JOBID
    exit
  endif
endif

if ($1 == 'remove') then  
  if (( $2 != 'R' ) && ( $2 != 'I' ) && ( $2 != 'H' )) then
    echo ""
    echo \'$2\' "is a wrong 'cwb_condor remove' option."
    echo "The availables options : cwb_condor remove I/R/H"
    echo ""
    exit
  else
    condor_q ${USER} -dag -nobatch -wide | grep -v "condor_dagman" | grep " $2 " | awk '{print $1}' | xargs condor_rm
    exit
  endif
endif

if ($1 == 'sdag') then  
  if ($2 == '') then
    echo ""
    echo "cwb_condor_sdag - Error : missing parameter"
    echo "Example : cwb_condor sdag nJOBS dag_file (optional)"
    echo "          where nJOBS is the number of jobs to be executed sequentialy on the nodes (>1)"
    echo ""
    exit
  else
    setenv CWB_CONDOR_NJOBS $2
    if ( $3 != '' ) then
      setenv CWB_CONDOR_DAG $3
    else
      unsetenv CWB_CONDOR_DAG
    endif
    root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_condor_sdag.C
    unsetenv CWB_CONDOR_NJOBS
    unsetenv CWB_CONDOR_DAG
    exit
  endif
endif

if ($1 == 'xdag') then  
  if ($2 == '' || $3 =='') then
    echo ""
    echo "cwb_condor_xdag - Error : missing parameter"
    echo "Example : cwb_condor xdag FRAME_GPS_START FRAME_LENGHT"
    echo "          See $CWB_MACROS/cwb_condor_xdag.C for more infos"
    echo ""
    exit
  else
    setenv CWB_CONDOR_FSTART $2
    setenv CWB_CONDOR_FLEN $3
    root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_condor_xdag.C
    unsetenv CWB_CONDOR_FSTART
    unsetenv CWB_CONDOR_FLEN
    exit
  endif
endif

if (($1 == 'create') || ($1 == 'recovery') || ($1 == 'resume')) then  
  if (($2 == '') || ($2 == 'FULL') || ($2 == 'INIT') || ($2 == 'STRAIN') || ($2 == 'CSTRAIN') || ($2 == 'COHERENCE') || ($2 == 'SUPERCLUSTER') || ($2 == 'LIKELIHOOD')) then   # $2 is a stage
    if ( $2 != '' ) then
      setenv CWB_STAGE_NAME $2
      if ( $3 != '' ) then		# define input directory
        setenv CWB_STAGE_INPUT $3
      else
        setenv CWB_STAGE_INPUT ""	# default output_dir is used as input dir
      endif
    else
      setenv CWB_STAGE_NAME "FULL"
    endif
    setenv CWB_STAGE_NAME "CWB_STAGE_"$CWB_STAGE_NAME
  else 
    echo "option : "$2" not allowed"
    echo "type cwb_condor to list the available options"
    exit
  endif
  # recovery but is done only if previous cwb_stage is present in the output dir
  if ($1 == 'resume') then  
    setenv CWB_STAGE_RESUME "TRUE"
  else 
    setenv CWB_STAGE_RESUME "FALSE"
  endif
endif

if ($1 == 'benchmark') then  
  if ($# > 1) then
    setenv CWB_BENCH_OPTS "$2"

    root -n -l ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_condor_benchmark.C

    unsetenv CWB_BENCH_OPTS
  else 
    $CWB_SCRIPTS/cwb_help.csh cwb_condor benchmark
    exit
  endif
endif

if ( $1 == 'create' ) then

  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_condor_create.C

endif

if ( $1 == 'submit' ) then

  if ( $2 != '' ) then
    setenv CWB_CONDOR_DAG $2
  else
    unsetenv CWB_CONDOR_DAG
  endif

  if (! $?_USE_PEGASUS ) then
    echo ""
  else
    voms-proxy-info --timeleft | awk '{if($1 < 600) exit 1}'	#_SKIP_CSH2SH_
    if ( $? != 0) then
      echo ""
      echo "voms-proxy timeleft too short, renew your proxy (Ex: voms-proxy-init -voms virgo)"
      echo ""
      voms-proxy-info --timeleft | awk 'BEGIN { OFS = ""; ORS = "" } ;{print "timeleft : "; print $1; print " sec\n"}'	#_SKIP_CSH2SH_
      echo ""
      exit
    endif

    $CWB_SCRIPTS/cwb_pegasus.csh status | tail -n 2 | head -n 1 | awk '{if($9 == "Running") exit 1}'	#_SKIP_CSH2SH_
    if ( $? != 0) then
      $CWB_SCRIPTS/cwb_pegasus.csh status	#_SKIP_CSH2SH_
      echo ""
      echo "jobs process already running ..."
      echo "to submit new jobs the current jobs must be removed : use cwb_premove"
      echo ""
      exit
    endif

    if ( $3 != '' ) then
      setenv CWB_PEGASUS_USITE $3
    else
      unsetenv CWB_PEGASUS_USITE 
    endif
  endif

  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_condor_submit.C

endif

if (($1 == 'resume') || ($1 == 'recovery')) then  

  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_condor_recovery.C

endif

if ( $1 == 'status' ) then
  $CWB_SCRIPTS/cwb_condor.csh status	#_SKIP_CSH2SH_
endif

if ( $1 == 'check' ) then

  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_condor_check.C

endif

if ( $1 == 'list' ) then

  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_condor_list.C

endif

if ( $1 == 'cleanup' ) then

  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_condor_cleanup.C

endif

# create cWB_analysis.log file
make -f $CWB_SCRIPTS/Makefile.log CMD_LINE="$cmd_line" svn >& /dev/null

unsetenv CWB_BENCHMARK_NAME
unsetenv CWB_STAGE_NAME
unsetenv CWB_STAGE_INPUT
unsetenv CWB_STAGE_RESUME
unsetenv CWB_CONDOR_DAG
unsetenv CWB_PEGASUS_USITE 

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

