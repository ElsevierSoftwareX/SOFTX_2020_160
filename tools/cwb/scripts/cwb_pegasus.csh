#!/bin/tcsh -f

if ($1 == '') then
  #$CWB_SCRIPTS/cwb_help.csh cwb_pegasus
  echo "available options :"
  echo ""
  echo "create jobs stuff      : create"
  echo "submit jobs            : submit"
  echo "recovery jobs          : recovery"
  echo "produce benchmark      : benchmark"
  echo "show status            : status"
  echo "re-submit failed jobs  : run"
  echo "delete all jobs        : remove"
  echo "show workflow analysis : analyzer"
  echo "make plots             : plots"
  echo "voms proxy timeleft    : timeleft"
  echo "voms-proxy-init        : voms"
  echo ""
  exit
endif

if (! $?_USE_PEGASUS ) then
  echo "environment _USE_PEGASUS not set"
  echo ""
  exit 
endif

if (( $1 != 'create' ) && ( $1 != 'submit' ) && ( $1 != 'recovery' ) && ( $1 != 'benchmark' ) && ( $1 != 'status' ) && ( $1 != 'run' ) && ( $1 != 'remove' ) && ( $1 != 'analyzer' ) && ( $1 != 'plots' ) && ( $1 != 'timeleft' ) && ( $1 != 'voms' )) then
  echo ""
  echo \'$1\' "is a wrong cwb_pegasus option"
  echo "type cwb_pegasus & return to see the availables options"
  echo ""
  exit
endif

setenv CWB_PEGASUS_OPTION $1 

if ( $CWB_PEGASUS_OPTION == 'voms' ) then
  if ( $2 != '' ) then
    setenv SITE_NAME $2
  else
    setenv SITE_NAME "virgo"
  endif

  if ( $3 != '' ) then
    setenv SITE_CATALOG $3
  else
    setenv SITE_CATALOG "virgo-voms.cnaf.infn.it"
  endif

  voms-proxy-init -voms $SITE_NAME -vomses $SITE_CATALOG
  exit
endif

if (($CWB_PEGASUS_OPTION == 'create') || ($CWB_PEGASUS_OPTION == 'recovery') || ($CWB_PEGASUS_OPTION == 'submit') || ($CWB_PEGASUS_OPTION == 'benchmarl')) then
  $CWB_SCRIPTS/cwb_condor.csh $CWB_PEGASUS_OPTION $2 $3
  exit
endif

voms-proxy-info --timeleft | awk '{if($1 <= 0) exit 1}'
if ( $? != 0) then
  echo ""
  echo "voms-proxy no timeleft, renew your proxy (Ex: voms-proxy-init -voms virgo)"
  echo ""
  exit
endif

if ( $CWB_PEGASUS_OPTION == 'timeleft' ) then
  voms-proxy-info --timeleft | awk 'BEGIN { OFS = ""; ORS = "" } ;{print "timeleft : "; print $1; print " sec\n"}'  
  exit
endif

if ( $CWB_PEGASUS_OPTION == 'status' ) then
  basename $PWD | awk 'BEGIN { OFS = ""; ORS = "" } ; {print "condor/"} ; {print $1} ' | xargs pegasus-status -l
  exit
endif

if ( $CWB_PEGASUS_OPTION == 'run' ) then
  basename $PWD | awk 'BEGIN { OFS = ""; ORS = "" } ; {print "condor/"} ; {print $1} ' | xargs pegasus-run
  exit
endif

if ( $CWB_PEGASUS_OPTION == 'remove' ) then
  basename $PWD | awk 'BEGIN { OFS = ""; ORS = "" } ; {print "condor/"} ; {print $1} ' | xargs pegasus-remove
  exit
endif

if ( $CWB_PEGASUS_OPTION == 'analyzer' ) then
  basename $PWD | awk 'BEGIN { OFS = ""; ORS = "" } ; {print "condor/"} ; {print $1} ' | xargs pegasus-analyzer
  exit
endif

if ( $CWB_PEGASUS_OPTION == 'plots' ) then
  basename $PWD | awk 'BEGIN { OFS = ""; ORS = "" } ; {print "condor/"} ; {print $1} ' | xargs pegasus-plots -p all
  exit
endif

exit
endif

