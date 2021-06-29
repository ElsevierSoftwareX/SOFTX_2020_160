#!/bin/tcsh -f

if ( $1 == '' ) then
  echo "specify command"
  echo "cwb_online option [create/file] python_file"
  echo "option=file: copy standard configuration in filename"
  echo "option=create: create directory from filename"
  exit
endif

#if ( ($1 != "create") && ($1 != "file") ) then
if ( ($1 != "create") ) then
  echo "cwb_online option [create/file] python_file"
  exit
endif

if  ( $1 == "create" ) then
  mkdir ${HOME}/tmp_ONLINE
  cp $2 ${HOME}/tmp_ONLINE/cWB_conf.py
  cp ${CWB_MACROS}/cwb_online.py ${HOME}/tmp_ONLINE/.
  cd ${HOME}/tmp_ONLINE
  ./cwb_online.py
  cd -
  rm -rf ${HOME}/tmp_ONLINE
endif

#if ( $1 == "file" ) then
#  cp ${CWB_ONLINE}/cWB_conf.py $2
#endif
