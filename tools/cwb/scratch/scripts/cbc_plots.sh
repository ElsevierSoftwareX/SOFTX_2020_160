if [[  -z $1  ]]; then
  echo ""
  echo 'cbc_plots mlabel bkg_tree(optional) ' 
  echo ""
  echo "            this command produces cbc plots + html page"
  echo "            when given the bkg tree it produces also FAD plots"
  echo ""
  echo "            mlabel : merge label"
  echo "    bkg_tree: entire path to the the tree root file   "
   echo '            Ex :   cbc_plots M1 ../TEST1/merge/wave_TEST1.M2.root '
  echo "  "
  return 0

else
  export CWB_MERGE_LABEL="$1"
 
  if [[  -n $2  ]]; then

  export CWB_BKG_TREE="$2"
  fi
fi

  #root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/../postproduction/cbc/cbc_plots.C
  root -n -l ${CWB_PPARMS_FILES} ${CWB_MACROS}/../postproduction/cbc/cbc_plots.C
  if [[  $? != 0 ]]; then
    echo ""
    echo "cbc_plots.C error : process terminated"
    echo ""
    return 0
   fi

  root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_mkhtml_index.C
  if [[  $? != 0 ]]; then
    echo ""
    echo "cwb_mkhtml_index.C error : process terminated"
    echo ""
    return 0
  fi

  root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/cwb_mkhtml_header.C
  if [[  $? != 0 ]]; then
    echo ""
    echo "cwb_mkhtml_header.C error : process terminated"
    echo ""
    return 0
  fi
  
  root -n -l -b ${CWB_PPARMS_FILES} ${CWB_MACROS}/../postproduction/cbc/cwb_mkhtml_cbc.C
  if [[  $? != 0 ]]; then
    echo ""
    echo "cwb_mkhtml_ebbh.C error : process terminated"
    echo ""
    return 0
  fi

  unset CWB_MERGE_LABEL
  unset CWB_BKG_TREE



  
