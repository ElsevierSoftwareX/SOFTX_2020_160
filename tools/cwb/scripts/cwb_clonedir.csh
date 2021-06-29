#!/bin/tcsh -f

onintr irq_ctrlc

set cmd_line="$0 $argv"

if (($1 == '') && ($2 == '')) then
  $CWB_SCRIPTS/cwb_help.csh cwb_clonedir
  exit
endif

if (($1 != '') && ($2 == '')) then
  echo ""
  echo "wrong cwb_clonedir parameters"
  echo "type cwb_clonedir & return to see the availables options"
  echo ""
  exit
endif

if ( ! -e $1 ) then 
  echo ""
  echo -n $1
  echo " not exist!"
  echo ""
  exit 
endif 

setenv CWB_CLONEDIR_PWD ${PWD}
setenv CWB_CLONEDIR_SRC $1
setenv CWB_CLONEDIR_OPTIONS "$3"

# extract dest dir
if ($2 == '.') then
  set dest_dir = `basename $CWB_CLONEDIR_SRC`
  setenv CWB_CLONEDIR_DEST $dest_dir
else
  setenv CWB_CLONEDIR_DEST $2
endif

if (( -e $CWB_CLONEDIR_DEST ) && ( "$3" == '')) then 
  echo ""
  echo -n $CWB_CLONEDIR_DEST
  echo " already exist!"
  echo ""
  exit 
endif 

# check options
if ( "$CWB_CLONEDIR_OPTIONS" != '' ) then
  setenv CWB_CLONEDIR_CHECK "CHECK"
  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_clonedir.C
  set return = $?
  if ( $return == 1 ) then
    echo "cwb_clonedir options are not correct"
    echo "type cwb_clonedir to see the available options"
    echo ""
    unsetenv CWB_CLONEDIR_CHECK
    unsetenv CWB_CLONEDIR_OPTIONS
    exit 
  endif
  unsetenv CWB_CLONEDIR_CHECK
endif

echo $CWB_CLONEDIR_SRC $CWB_CLONEDIR_DEST
if ( -e $CWB_CLONEDIR_DEST ) then
  # execute options
  cd $CWB_CLONEDIR_DEST
  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_clonedir.C
  if ( $? != 0 ) then
    echo ""
    echo "Error executing cwb_clonedir"
    echo ""
  endif 
  cd $CWB_CLONEDIR_PWD
  unsetenv CWB_CLONEDIR_OPTIONS
  unsetenv CWB_CLONEDIR_PWD
  unsetenv CWB_CLONEDIR_SRC
  unsetenv CWB_CLONEDIR_DEST
  exit
else
  # create dest dir
  mkdir $CWB_CLONEDIR_DEST
endif

# copy src dir to dest dir
if ( -e $CWB_CLONEDIR_SRC/README ) then
  cp -r $CWB_CLONEDIR_SRC/README $CWB_CLONEDIR_DEST/.
endif
if ( -e $CWB_CLONEDIR_SRC/config ) then
  cp -r $CWB_CLONEDIR_SRC/config $CWB_CLONEDIR_DEST/.
  # makes ROOT5 macro configuration files ROOT6 compliant 
  root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_MACROS}/FixUserParametersROOT6.C\(\"$CWB_CLONEDIR_DEST/config/user_parameters.C\",\"p\"\)
  root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_MACROS}/FixUserParametersROOT6.C\(\"$CWB_CLONEDIR_DEST/config/user_pparameters.C\",\"pp\"\)
endif
if ( -e $CWB_CLONEDIR_SRC/input ) then
  cp -r $CWB_CLONEDIR_SRC/input  $CWB_CLONEDIR_DEST/.
endif
if ( -e $CWB_CLONEDIR_SRC/macro ) then
  cp -r $CWB_CLONEDIR_SRC/macro  $CWB_CLONEDIR_DEST/.
  # force plugin to be recompiled
  touch $CWB_CLONEDIR_DEST/macro/*_C.so -d "10 years ago"
endif

# create full dest dir 
cd $CWB_CLONEDIR_DEST
setenv CWB_MKDIR_WRKDIR ""
setenv CWB_MKDIR_OPTION "batch"
root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_mkdir.C
if ( $? != 0) then
  echo ""
  cd $CWB_CLONEDIR_PWD
  echo "cwb_clonedir.C error : process terminated"
  echo ""
  echo "remove directory : " $CWB_CLONEDIR_DEST 
  echo ""
  exit
endif
#root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_condor_create.C

# apply options
if ( "$CWB_CLONEDIR_OPTIONS" != '' ) then
  root -n -l -b ${CWB_PARMS_FILES} ${CWB_MACROS}/cwb_clonedir.C
endif

# create cWB_analysis.log file
make -f $CWB_SCRIPTS/Makefile.log CMD_LINE="$cmd_line" svn >& /dev/null

cd $CWB_CLONEDIR_PWD

echo ""
echo "The new working dir is : " $CWB_CLONEDIR_DEST
echo ""

# remove envs
unsetenv CWB_CLONEDIR_OPTIONS
unsetenv CWB_CLONEDIR_PWD
unsetenv CWB_CLONEDIR_SRC
unsetenv CWB_CLONEDIR_DEST

exit 0
irq_ctrlc:
  ps T | grep root | awk '{print $1}' | xargs kill -9
  exit 1

