#!/bin/tcsh -f

if ($1 == '') then
  #$CWB_SCRIPTS/cwb_help.csh cwb_lsf
  echo "available options :"
  echo ""
  echo "create jobs stuff      : create"
  echo "submit jobs            : submit"
  echo "recovery jobs          : recovery"
  echo "produce benchmark      : benchmark"
  echo "show queue status      : queue"
  echo "show all jobs status   : status"
  echo "show job log           : status jobID"
  echo "dump job log           : log jobName"
  echo "delete all jobs        : kill"
  echo "delete job             : kill jobID"
  echo "suspend all jobs       : stop"
  echo "suspend job            : stop jobID"
  echo "resume all jobs        : resume"
  echo "resume job             : resume jobID"
  echo "uncompress tgz files   : utar dir_name"
  echo ""
  exit
endif

if (! $?_USE_LSF ) then
  echo "environment _USE_LSF not set"
  echo ""
  exit 
endif

if (( $1 != 'create' ) && ( $1 != 'submit' ) && ( $1 != 'recovery' ) && ( $1 != 'benchmark' ) && ( $1 != 'queue' ) && ( $1 != 'status' ) && ( $1 != 'resume' ) && ( $1 != 'kill' ) && ( $1 != 'stop' ) && ( $1 != 'resume' ) && ( $1 != 'log' ) && ( $1 != 'utar' ) ) then
  echo ""
  echo \'$1\' "is a wrong cwb_lsf option"
  echo "type cwb_lsf & return to see the availables options"
  echo ""
  exit
endif

setenv CWB_LSF_OPTION $1 

if (($CWB_LSF_OPTION == 'create') || ($CWB_LSF_OPTION == 'recovery') || ($CWB_LSF_OPTION == 'submit') || ($CWB_LSF_OPTION == 'benchmark') || ($CWB_LSF_OPTION == 'resume')) then
  $CWB_SCRIPTS/cwb_condor.csh $CWB_LSF_OPTION $2 $3
  exit
endif

if ( $CWB_LSF_OPTION == 'queue' ) then
  echo ""
  bqueues $LSF_QUEUE
  echo ""
  exit
endif

if ( $CWB_LSF_OPTION == 'status' ) then
  if ( $2 != '' ) then
    bpeek -f $2
    exit
  else 
    basename $PWD | awk -v user="$USER" 'BEGIN { OFS = ""; ORS = "" } ; {print "/"} ; {print user} ; {print "/"} ; {print $1} ' | xargs bjobs -g
    exit
  endif
endif

if ( $CWB_LSF_OPTION == 'resume' ) then
  if ( $2 != '' ) then
    bresume $2
    exit
  else 
    basename $PWD | awk -v user="$USER" 'BEGIN { OFS = ""; ORS = "" } ; {print "/"} ; {print user} ; {print "/"} ; {print $1} ; {print "0"} ' | xargs bresume -g
    exit
  endif
endif

if ( $CWB_LSF_OPTION == 'stop' ) then
  if ( $2 != '' ) then
    bstop $2
    exit
  else 
    basename $PWD | awk -v user="$USER" 'BEGIN { OFS = ""; ORS = "" } ; {print "/"} ; {print user} ; {print "/"} ; {print $1} ; {print "0"} ' | xargs bstop -g
    exit
  endif 
endif

if ( $CWB_LSF_OPTION == 'kill' ) then
  if ( $2 != '' ) then
    bkill -r $2
    exit
  else 
    basename $PWD | awk -v user="$USER" 'BEGIN { OFS = ""; ORS = "" } ; {print "/"} ; {print user} ; {print "/"} ; {print $1} ; {print " 0"} ' | xargs bkill -g
    exit
  endif
endif

if ( $CWB_LSF_OPTION == 'log' ) then
  #setenv CWB_LSF_JOBID $2 
  #basename $PWD | awk -v user="$USER" -v jobid="$CWB_LSF_JOBID" 'BEGIN { OFS = ""; ORS = "" } ; {print jobid} ; {print "_"} ; {print $1} ' | xargs bpeek -J 
  bpeek $2
  exit
endif

if ( $CWB_LSF_OPTION == 'utar' ) then
  if ( $2 != '' ) then
    setenv CWB_LSF_UTAR_DIR $2
    bash -c 'cd $CWB_LSF_UTAR_DIR; for a in `ls -1 *.tgz` ; do tar xzf $a ; rm $a; done'
  else
    echo "Error : the directory is not defined"
    echo ""
  endif
  unsetenv CWB_LSF_UTAR_DIR
  exit
endif

exit
