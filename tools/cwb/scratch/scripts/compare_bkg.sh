if [[  -z $1 || -z $2  ]]; then
  echo ""
  echo 'compare_bkg mlabel ntuple' 
  echo ""
  echo "            this command does some comparisons between the rho distributions of two productions "
  echo "            the rho distribution of mlabel X from the current directory"
  echo "            vs the corresponding one of the ntuple Y "
  echo "            mlabel : merge label"
  echo "            ntuple : compare ntuple"
  echo " "
  echo '            Ex :   compare_bkg M1 /home/salemi/NINJA2/test_bkg_2G_5/merge/test_bkg_2G_5.M2.root '
  echo "  "
  return 0
fi

export CWB_MERGE_LABEL="$1"
export CWB_COMPARE_TREE="$2"

root -n -b -l ${CWB_PPARMS_FILES} ${CWB_MACROS}/../postproduction/burst/compare_bkg.C
unset CWB_MERGE_LABEL
unset CWB_COMPARE_TREE
