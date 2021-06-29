#!/bin/bash -f
# this script is used to create the WWW directory structure for the cWB report access from www

if [[ -z $SITE_CLUSTER  ]]; then
  echo "cwb_create_www.sh error ..."
  echo "directory $SITE_CLUSTER is not defined !!!" 
  exit 1
fi

if [[ "$SITE_CLUSTER" != CIT && "$SITE_CLUSTER" != VIRTUALBOX && "$SITE_CLUSTER" != DOCKER && "$SITE_CLUSTER" != USER  ]]; then
  echo "cwb_create_www.sh error ..."
  echo "directory $SITE_CLUSTER is not allowed !!!" 
  exit 1
fi

if [[ "$SITE_CLUSTER" = VIRTUALBOX ]]; then

  CWB_WWW="/home/cwb/waveburst/WWW"
  CWB_LIB="/home/cwb/waveburst/git/cWB/library"
  CWB_DOC="/home/cwb/waveburst//git/cWB/documentation"
  ALADIN="/home/cwb/waveburst/SOFT/SCRIPTS/aladin-1.0"
  SHADOWBOX="/home/cwb/waveburst/SOFT/SCRIPTS/shadowbox-3.0.3"

fi

if [[ "$SITE_CLUSTER" = CIT ]]; then

  CWB_WWW="$HOME/public_html/waveburst/WWW"
  CWB_LIB="$HOME/git/cWB/library"
  CWB_DOC="$HOME/git/cWB/documentation"
  ALADIN="$HOME/SOFT/SCRIPTS/aladin-1.0"
  SHADOWBOX="$HOME/SOFT/SCRIPTS/shadowbox-3.0.3"

fi

if [[ "$SITE_CLUSTER" = DOCKER  || "$SITE_CLUSTER" = USER ]]; then

  CWB_WWW="$HOME/WWW"
  CWB_LIB="$HOME/git/cWB/library"
  CWB_DOC="$HOME/git/cWB/documentation"
  ALADIN="$HOME/SOFT/SCRIPTS/aladin-1.0"
  SHADOWBOX="$HOME/SOFT/SCRIPTS/shadowbox-3.0.3"

fi

if [[ -d $CWB_WWW ]]; then
  echo "cwb_create_www.sh error ... "
  echo "directory $CWB_WWW already exist !!!" 
  exit 1
fi

#create WWW dir
mkdir $CWB_WWW
cd $CWB_WWW

# create link to banner
ln -s $CWB_LIB/html banner

# create link to logo
ln -s $CWB_LIB/tools/cwb/www/logo logo

# make ced dir
mkdir ced
cd ced

# create link to ced index files
ln -s $CWB_LIB/tools/cwb/www/index index

# create ced script dir
mkdir scripts
cd scripts

# create link to ced scripts dir
ln -s $CWB_LIB/tools/cwb/www ced

# create link to aladin 
# not mandatory
if [[ -d $ALADIN ]]; then
  ln -s $ALADIN aladin
fi

# create link to shadowbox 
# not mandatory
if [[ -d $SHADOWBOX ]]; then
  ln -s $SHADOWBOX shadowbox
fi

# create root dir
cd ../..
mkdir root
cd root
# create link to root scripts
ln -s $CWB_LIB/html/etc/html scripts

#create doc dir
cd ..
mkdir doc
cd doc
# create cwb doc dir
mkdir cwb
cd cwb
# create link to cwb manual
ln -s $CWB_DOC/_build/html man
cd ../..

