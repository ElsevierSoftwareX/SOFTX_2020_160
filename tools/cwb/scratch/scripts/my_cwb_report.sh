if [[ -z $1 ]]; then
  echo ""
  echo "cwb_report label[create/publish/remove/clean]/list[report/merge] option Xlag Xslag (only for production)" 
  echo ""
  echo " default Xlag=-1    : the selected lags are  : lag>0"
  echo " default Xslag=-1   : the selected slags are : slag>=0"
  echo ""
  echo " Xlag=-1 & Xslag=-1 : lag>0 & slag>=0       (BACKGROUND)"
  echo " Xlag =0 & Xslag =0 : lag=0 & slag=0        (ZERO LAG !!!)"
  echo " Xlag>=0 & Xslag>=0 : lag=Xlag & slag=Xslag"
  echo " Xlag>=0 & Xslag=-1 : lag=Xlag & slag=any"
  echo " Xlag=-1 & Xslag>=0 : lag=any & slag=Xslag"
  echo ""
  return 0
fi

if [[ $1 = "list" ]]; then
  if [[ $2 = "merge" ]]; then
    root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE} ${CWB_MACROS}/cwb_dump_merge_dir.C
  else 
    root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE} ${CWB_MACROS}/cwb_dump_report_dir.C
  fi 
  return 0
fi

export CWB_MERGE_LABEL=$1

if [[ $2 = "remove" ]]; then
  root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE} ${CWB_MACROS}/cwb_rm_report_dir.C
  if [[  $? != 0 ]]; then
    echo ""
    echo "cwb_rm_report_dir.C error : process terminated"
    echo ""
    return 0 
  fi
  root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE} ${CWB_MACROS}/cwb_dump_report_dir.C
  return 0
fi

if [[  $2 = "publish"  ||  $2 = "clean"  ]]; then
  export CWB_PP_DATA_DIR=$1
  export CWB_PUBLISH_OPTION=$2
  root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE} ${CWB_MACROS}/cwb_publish.C
  if [[  $? != 0 ]]; then
    echo ""
    echo "cwb_publish.C error : process terminated"
    echo ""
    return 0 
  fi
  root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE} ${CWB_MACROS}/cwb_dump_report_dir.C
  return 0
fi

if [[ $2 != "create" ]]; then
  echo ""
  echo "cwb_report label/list[report/merge] create/publish/remove/clean  option Xlag(default lags>0 - prod)"
  echo ""
  return 0
fi

root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE} ${CWB_MACROS}/cwb_getsim.C

# simulation

if [[  $? = 1  ]]; then
  export CWB_LAG_NUMBER=-1
  export CWB_SLAG_NUMBER=-1

  root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE} ${CWB_PPARAMETERS_FILE} ${CWB_UPPARAMETERS_FILE} ${CWB_EPPARAMETERS_FILE} ${CWB_MACROS}/cwb_simplot.C
  if [[  $? != 0 ]]; then
    echo ""
    echo "cwb_simplot.C error : process terminated"
    echo ""
    return 0 
  fi
  root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE} ${CWB_PPARAMETERS_FILE} ${CWB_UPPARAMETERS_FILE} ${CWB_EPPARAMETERS_FILE} ${CWB_MACROS}/cwb_mkeff.C
  if [[  $? != 0 ]]; then
    echo ""
    echo "cwb_mkeff.C error : process terminated"
    echo ""
    return 0 
  fi
  root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE} ${CWB_PPARAMETERS_FILE} ${CWB_UPPARAMETERS_FILE} ${CWB_EPPARAMETERS_FILE} ${CWB_MACROS}/cwb_mkhtml_index.C
  if [[  $? != 0 ]]; then
    echo ""
    echo "cwb_mkhtml_index.C error : process terminated"
    echo ""
    return 0 
  fi
  root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE} ${CWB_PPARAMETERS_FILE} ${CWB_UPPARAMETERS_FILE} ${CWB_EPPARAMETERS_FILE} ${CWB_MACROS}/cwb_mkhtml_header.C
  if [[  $? != 0 ]]; then
    echo ""
    echo "cwb_mkhtml_header.C error : process terminated"
    echo ""
    return 0 
  fi
  root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE} ${CWB_PPARAMETERS_FILE} ${CWB_UPPARAMETERS_FILE} ${CWB_EPPARAMETERS_FILE} ${CWB_MACROS}/cwb_mkhtml_sim.C
  if [[  $? != 0 ]]; then
    echo "cwb_mkhtml_sim.C error : process terminated"
    return 0 
  fi

  return 0
fi

#production

if [[  $? = 0  ]]; then
  if [[ -z $3 ]]; then
    export CWB_LAG_NUMBER=-1
  else  
    export CWB_LAG_NUMBER=$3
  fi

  if [[ $4 = "" ]]; then
    export CWB_SLAG_NUMBER=-1
  else  
    export CWB_SLAG_NUMBER=$4
  fi

  root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE} ${CWB_PPARAMETERS_FILE} ${CWB_UPPARAMETERS_FILE} ${CWB_EPPARAMETERS_FILE} ${CWB_MACROS}/cwb_netplot_gen_1.C
  if [[  $? != 0 ]]; then
    echo ""
    echo "cwb_netplot_gen_1.C error : process terminated"
    echo ""
    return 0 
  fi
  root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE} ${CWB_PPARAMETERS_FILE} ${CWB_UPPARAMETERS_FILE} ${CWB_EPPARAMETERS_FILE} ${CWB_MACROS}/cwb_netplot_gen_2.C
  if [[  $? != 0 ]]; then
    echo "cwb_netplot_gen_2.C error : process terminated"
    return 0 
  fi
  root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE} ${CWB_PPARAMETERS_FILE} ${CWB_UPPARAMETERS_FILE} ${CWB_EPPARAMETERS_FILE} ${CWB_MACROS}/lag.C
  if [[  $? != 0 ]]; then
    echo ""
    echo "lag.C error : process terminated"
    echo ""
    return 0
  fi
  root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE} ${CWB_PPARAMETERS_FILE} ${CWB_UPPARAMETERS_FILE} ${CWB_EPPARAMETERS_FILE} ${CWB_MACROS}/slag.C
  if [[  $? != 0 ]]; then
    echo "slag.C error : process terminated"
    return 0
  fi
  root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE} ${CWB_PPARAMETERS_FILE} ${CWB_UPPARAMETERS_FILE} ${CWB_EPPARAMETERS_FILE} ${CWB_MACROS}/cwb_mkhtml_index.C
  if [[  $? != 0 ]]; then
    echo ""
    echo "cwb_mkhtml_index.C error : process terminated"
    echo ""
    return 0 
  fi
  root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE} ${CWB_PPARAMETERS_FILE} ${CWB_UPPARAMETERS_FILE} ${CWB_EPPARAMETERS_FILE} ${CWB_MACROS}/cwb_mkhtml_header.C
  if [[  $? != 0 ]]; then
    echo ""
    echo "cwb_mkhtml_header.C error : process terminated"
    echo ""
    return 0 
  fi
  root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE} ${CWB_PPARAMETERS_FILE} ${CWB_UPPARAMETERS_FILE} ${CWB_EPPARAMETERS_FILE} ${CWB_MACROS}/cwb_mkhtml_prod.C
  if [[  $? != 0 ]]; then
    echo ""
    echo "cwb_mkhtml_prod.C error : process terminated"
    echo ""
    return 0 
  fi
  root -n -l -b ${CWB_ROOTLOGON_FILE} ${CWB_PARAMETERS_FILE} ${CWB_UPARAMETERS_FILE} ${CWB_EPARAMETERS_FILE} ${CWB_PPARAMETERS_FILE} ${CWB_UPPARAMETERS_FILE} ${CWB_EPPARAMETERS_FILE} ${CWB_MACROS}/cwb_condor_create_ced.C
  if [[  $? != 0 ]]; then
    echo ""
    echo "cwb_condor_create_ced.C error : process terminated"
    echo ""
    return 0
  fi

  return 0
fi

