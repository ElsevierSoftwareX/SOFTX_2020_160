#!/bin/bash

if [[ -z $1 ]]; then
  echo ""
  echo "cwb_pegasus_create.sh :  missing dag file"
  echo ""
  exit 1
else 
  if [ -f $1 ];
  then
    CWB_DAGFILE=$1                              # DAG FILE
    echo "input  dag file : " \'$CWB_DAGFILE\' 
    CWB_DAXFILE="${CWB_DAGFILE/.dag/.dax}"      # DAX FILE
    echo "output dax file : " \'$CWB_DAXFILE\' 
    CWB_OUTFILE="${CWB_DAGFILE/.dag/.out}"      # OUT FILE
    echo "output out file : " \'$CWB_OUTFILE\' 
    CWB_INFILE="${CWB_DAGFILE/.dag/.in}"        # IN  FILE : in tgz list of files 
    echo "output in file  : " \'$CWB_INFILE\' 
    CWB_TGZFILE="${CWB_DAGFILE/.dag/.tgz}"      # TGZ FILE  
    echo "output tgz file : " \'$CWB_TGZFILE\' 
    if [ -f $CWB_DAXFILE ];
    then
      echo "output dax file : " \'$CWB_DAXFILE\' "already exist !!!"
      read -p "do you want to overwrite it ? (y/n) " 
      if [ "$REPLY" != "y" ];
      then       
        exit 0
      fi
    fi
  else
    echo "File $1 does not exist."
    exit 1
  fi
fi

set -e

tar czvfh $CWB_TGZFILE -T $CWB_INFILE 2>&1 >/dev/null

cd ..
TOPDIR=`pwd`
TOPDIRNAME=${PWD##*/}
cd -

# pegasus bin directory is needed to find keg
BIN_DIR=`pegasus-config --bin`

# build the dax generator
CLASSPATH=`pegasus-config --classpath`
export CLASSPATH=".:$CLASSPATH"
javac ${CWB_SCRIPTS}/cwb_dag2dax.java -d .

# generate the dax
echo "output watenv file : " \'$CWB_PEGASUS_WATENV\' 
java cwb_dag2dax $BIN_DIR $CWB_DAGFILE $CWB_TGZFILE $CWB_DAXFILE $CWB_OUTFILE $CWB_PEGASUS_WATENV

# create the site catalog
sed "s|WORKDIR|$TOPDIR|g" $HOME_WAT/tools/pegasus.sites > sites.xml3

