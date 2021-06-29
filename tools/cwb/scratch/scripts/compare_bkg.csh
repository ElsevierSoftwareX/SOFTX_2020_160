#!/bin/tcsh -f

if ( $1 == '' || $2 == '' ) then
  echo ""
  echo 'compare_bkg mlabel ntuple' 
  echo ""
  echo "            this command does some comparisons between the rho distributions of two productions "
  echo "            the rho distribution of mlabel X from the current directory"
  echo "            vs the corresponding one of the ntuple Y	"
  echo "            mlabel : merge label"
  echo "            ntuple : compare ntuple"
  echo " "
  echo '            Ex :   compare_bkg M1 /home/salemi/NINJA2/test_bkg_2G_5/merge/test_bkg_2G_5.M2.root '
  echo "  "
  exit
endif

setenv CWB_MERGE_LABEL "$1"
setenv CWB_COMPARE_TREE "$2"

 root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/../postproduction/burst/compare_bkg.C

 root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_mkhtml_index.C
  if ( $? != 0) then
    echo ""
    echo "cwb_mkhtml_index.C error : process terminated"
    echo ""
    exit
  endif

 

 # root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_mkhtml_header.C
  if ( $? != 0) then
    echo ""
    echo "cwb_mkhtml_header.C error : process terminated"
    echo ""
    exit
  endif
  
 # root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/../postproduction/cbc/cwb_mkhtml_cbc.C
  if ( $? != 0) then
    echo ""
    echo "cwb_mkhtml_ebbh.C error : process terminated"
    echo ""
    exit
  endif


