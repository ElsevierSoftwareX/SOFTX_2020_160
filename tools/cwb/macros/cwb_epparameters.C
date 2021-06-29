// exec pparameters macro, used after the user_pparameters.C

{
  CheckAnalysis();         // check consistency with the environment CWB_ANALYSIS

  if(simulation<0) return; // super report mode

  // ----------------------------------------------------------
  // Warning if zero lag analysis
  // ----------------------------------------------------------

  if(!pp_batch&&(cwb_lag_number==0)&&(cwb_slag_number==0)&&(simulation==0)) {
    char answer[1024];
    strcpy(answer,"");
    do {
      cout << "You are going to produce the report for lag zero" << endl;
      cout << "Do you want to procede  ? (y/n) ";
      cin >> answer;
      cout << endl << endl;
    } while ((strcmp(answer,"y")!=0)&&(strcmp(answer,"n")!=0));
    if (strcmp(answer,"n")==0) exit(1);
  }

  // ----------------------------------------------------------
  // Declare standard cuts & standard pp_label
  // ----------------------------------------------------------

  #include "GToolbox.hh"

  if(pp_rho_min<0.) {cout << "Error : pp_rho_min must be >=0" << endl;exit(1);}
  if(pp_rho_max<=pp_rho_min) {cout << "Error : pp_rho_max must be > pp_rho_min" << endl;exit(1);}

  pp_rho_bin = (pp_rho_max-pp_rho_min)/pp_drho;

  if(pp_rho_bin>100) pp_rho_bin=100;	// the upper value is limited because it requires big memory resources 
  if(pp_rho_bin<=0) {cout << "Error : pp_rho_bin must be >0, check pp_rho_min, pp_rho_max, pp_drho" << endl;exit(1);}

  if(TString(gSystem->Getenv("CWB_REPORT_PE"))!="") {		// parameter estimation report
    if(optim) pp_label=TString("_PE_")+TString(SEARCH())+"SRA";
    else      pp_label=TString("_PE_")+TString(SEARCH())+"MRA";
  } else if(TString(gSystem->Getenv("CWB_REPORT_CBC"))!="") {	// CBC report
    if(optim) pp_label=TString("_CBC_")+TString(SEARCH())+"SRA";
    else      pp_label=TString("_CBC_")+TString(SEARCH())+"MRA";
  } else {							// standard report
    if(optim) pp_label=TString("_")+TString(SEARCH())+"SRA";	
    else      pp_label=TString("_")+TString(SEARCH())+"MRA";
  }

  double ifLow=fLow;
  double ifHigh=fHigh;

#ifdef PP_64_200_HZ
  //pp_label="64_200";
  ifLow=64;ifHigh=200; 
#endif
#ifdef PP_200_2048_HZ
  //pp_label="200_2048";
  ifLow=200;ifHigh=2048; 
#endif
#ifdef PP_1600_5000_HZ
  //pp_label="1600_5000";
  ifLow=1600;ifHigh=5000;
#endif

#ifdef CAT2_VETO
  pp_label=pp_label+"_cat2";
#endif
#ifdef HVETO_VETO
  pp_label=pp_label+"_hveto";
#endif
#ifdef CAT3_VETO
  pp_label=pp_label+"_cat3";
#endif
#ifdef PEM_VETO
  pp_label=pp_label+"_pem";
#endif

  char L_cor[32]=""; 
  char L_cut[32]=""; 
  char L_CUT[32]=""; 
  char L_vED[32]=""; 
  char L_scc[32]=""; 
  char L_pen[32]=""; 
  char L_win[32]=""; 
  char L_ifar[32]=""; 

  if(T_cor>=0) sprintf(L_cor,"i%dcc%.2i", pp_inetcc, (int)(T_cor*100));
  if(T_cut>=0) sprintf(L_cut,"i%drho%i", pp_irho, (int)(T_cut*10));   
  if(T_vED>0)  sprintf(L_vED,"vED%.2i",(int)(T_vED*100));  
  if(T_scc>0)  sprintf(L_scc,"scc%.2i",(int)(T_scc*100)); 
  if(T_pen>0)  sprintf(L_pen,"pen%.2i",(int)(T_pen*100)); 
  if(simulation) {if(T_win>0) sprintf(L_win,"win%.2i",(int)(T_win*10)); 
                  else        {cout << "Error : T_win must be >0" << endl;exit(1);}}
  TString sT_ifar = TString::Format("%g",T_ifar);
  sT_ifar.ReplaceAll(".","d");
  if(T_ifar>0) sprintf(L_ifar,"ifar%s",sT_ifar.Data()); 

  if(TString(L_cor)!="") pp_label=pp_label+TString("_"+TString(L_cor)); 
  if(TString(L_cut)!="") pp_label=pp_label+TString("_"+TString(L_cut)); 
  if(TString(L_CUT)!="") pp_label=pp_label+TString("_"+TString(L_CUT)); 
  if(TString(L_vED)!="") pp_label=pp_label+TString("_"+TString(L_vED)); 
  if(TString(L_scc)!="") pp_label=pp_label+TString("_"+TString(L_scc)); 
  if(TString(L_pen)!="") pp_label=pp_label+TString("_"+TString(L_pen)); 
  if(TString(L_win)!="") pp_label=pp_label+TString("_"+TString(L_win)); 
  if(TString(L_ifar)!="") pp_label=pp_label+TString("_"+TString(L_ifar)); 

  if(user_pp_label.Sizeof()>1) pp_label=pp_label+TString("_")+user_pp_label;

  // ----------------------------------------------------------
  // check pp parameters
  // ----------------------------------------------------------

  if((pp_irho!=0)&&(pp_irho!=1)) {
    cout << "cwb_epparameters.C : Error - wrong pp_irho=" << pp_irho << " value "
         << "must be 0/1" << endl;
    gSystem->Exit(1);  
  }
  if((pp_inetcc!=0)&&(pp_inetcc!=1)) {
    cout << "cwb_epparameters.C : Error - wrong pp_irho=" << pp_inetcc << " value "
         << "must be 0/1" << endl;
    gSystem->Exit(1);  
  }

  // ----------------------------------------------------------
  // frequency cuts
  // ----------------------------------------------------------

#ifdef PP_64_200_HZ
  lowFCUT[nFCUT] = 200.;  highFCUT[nFCUT] = 2048.; nFCUT++;
#endif
#ifdef PP_200_2048_HZ
  lowFCUT[nFCUT] = 0.;    highFCUT[nFCUT] = 200.;  nFCUT++;
#endif
#ifdef PP_1600_5000_HZ
  lowFCUT[nFCUT] = 0.;    highFCUT[nFCUT] = 1600.; nFCUT++;
  lowFCUT[nFCUT] = 5000.; highFCUT[nFCUT] = 8192.; nFCUT++;
#endif

  // check frequency cuts consistency
  for(int i=0;i<nFCUT;i++) {
    if(highFCUT[i]<lowFCUT[i]) {
      cout<< "Error in Frequency cuts : lowFCUT imust be < highFCUT" 
          << lowFCUT[i] << " " << highFCUT[i] << endl;
      exit(1);
    }
  }
  // select freequncy range for plots
  double freqLow=fLow;
  double freqHigh=fHigh;
  for(int i=0;i<nFCUT;i++) {
    if(fLow>=lowFCUT[i] && fLow<=highFCUT[i])   freqLow=highFCUT[i];
    if(fHigh>=lowFCUT[i] && fHigh<=highFCUT[i]) freqHigh=lowFCUT[i];
  }
#ifdef _USE_ROOT6
  // export to CINT (this is a fix for ROOT6)
  char epp_cmd[128]; 
  sprintf(epp_cmd,"freqLow  = %f;",freqLow);  EXPORT(double,freqLow,epp_cmd)
  sprintf(epp_cmd,"freqHigh = %f;",freqHigh); EXPORT(double,freqHigh,epp_cmd)
#endif

  // union of freq rages to get label
  waveSegment seg;
  vector<waveSegment> iseg; iseg.clear();
  seg.index=0; seg.start=0; seg.stop=fLow; iseg.push_back(seg);
  seg.index=1; seg.start=fHigh; seg.stop=inRate/2.; iseg.push_back(seg);
  for(int i=0;i<nFCUT;i++) {
    seg.index=i+2; seg.start=lowFCUT[i]; seg.stop=highFCUT[i]; 
    iseg.push_back(seg);
  }
  vector<waveSegment> oseg = CWB::Toolbox::unionSegments(iseg);
  char freqs_label[1024]="";
  for(int n=0;n<oseg.size()-1;n++) {
    sprintf(freqs_label,"%s_freq%d_%d",freqs_label,int(oseg[n].stop),int(oseg[n+1].start));
  }
  TString sfreqs_label = freqs_label;
  sfreqs_label.ReplaceAll(".","d");
  pp_label = pp_label+sfreqs_label;

  // ----------------------------------------------------------
  // Build title & subtitle for reports
  // ----------------------------------------------------------

  strcpy(RunLabel,RUN_LABEL);

  char xtitle[1024];
  if(strlen(ifo[0])>0) strcpy(xtitle,ifo[0]); else strcpy(xtitle,detParms[0].name);
  for(int i=1;i<nIFO;i++) {
    TString xtmp = xtitle;
    if(strlen(ifo[i])>0) sprintf(xtitle,"%s-%s",xtmp.Data(),ifo[i]);           // built in detector
    else                 sprintf(xtitle,"%s-%s",xtmp.Data(),detParms[i].name); // user define detector
  }

  if(simulation) {
    sprintf(title,"SIMULATION - Network %s [%4.1f-%4.1f]",xtitle,ifLow,ifHigh);
  } else {
    if((cwb_lag_number==-1)&&(cwb_slag_number==-1)) {
      sprintf(title,"BACKGROUND NON ZERO LAG  - Network %s [%4.1f-%4.1f]",xtitle,ifLow,ifHigh);
    } 
    if((cwb_lag_number>=0)&&(cwb_slag_number>=0)) {
      sprintf(title,"BACKGROUND LAG %d SLAG %d - Network %s [%4.1f-%4.1f]",cwb_lag_number,cwb_slag_number,xtitle,ifLow,ifHigh);
    }
    if((cwb_lag_number>=0)&&(cwb_slag_number==-1)) {
      sprintf(title,"BACKGROUND LAG %d - Network %s [%4.1f-%4.1f]",cwb_lag_number,xtitle,ifLow,ifHigh);
    }
    if((cwb_lag_number==-1)&&(cwb_slag_number>=0)) {
      sprintf(title,"BACKGROUND SLAG %d - Network %s [%4.1f-%4.1f]",cwb_slag_number,xtitle,ifLow,ifHigh);
    }
  }


  // ----------------------------------------------------------
  // Build output directory names
  // ----------------------------------------------------------

  if(cwb_merge_label.Sizeof()>1) sprintf(pp_dir,"%s/%s",pp_dir,cwb_merge_label.Data());
  if(pp_label.Sizeof()>1) {
    if((cwb_lag_number==-1)&&(cwb_slag_number==-1)) {
      sprintf(pp_dir,"%s.R%s",pp_dir,pp_label.Data());
    } 
    if((cwb_lag_number>=0)&&(cwb_slag_number>=0)) {
      sprintf(pp_dir,"%s.R%s_lag%d_slag%d",pp_dir,pp_label.Data(),cwb_lag_number,cwb_slag_number);
    }
    if((cwb_lag_number>=0)&&(cwb_slag_number==-1)) {
      sprintf(pp_dir,"%s.R%s_lag%d",pp_dir,pp_label.Data(),cwb_lag_number);
    }
    if((cwb_lag_number==-1)&&(cwb_slag_number>=0)) {
      sprintf(pp_dir,"%s.R%s_slag%d",pp_dir,pp_label.Data(),cwb_slag_number);
    }
  } else {
    if((cwb_lag_number==-1)&&(cwb_slag_number==-1)) {
      sprintf(pp_dir,"%s.R",pp_dir);
    } 
    if((cwb_lag_number>=0)&&(cwb_slag_number>=0)) {
      sprintf(pp_dir,"%s.Rlag%d_slag%d",pp_dir,cwb_lag_number,cwb_slag_number);
    }
    if((cwb_lag_number>=0)&&(cwb_slag_number==-1)) {
      sprintf(pp_dir,"%s.Rlag%d",pp_dir,cwb_lag_number);
    }
    if((cwb_lag_number==-1)&&(cwb_slag_number>=0)) {
      sprintf(pp_dir,"%s.Rslag%d",pp_dir,cwb_slag_number);
    }
  }
  sprintf(netdir,"%s/%s",pp_dir,pp_data_dir);

  // ----------------------------------------------------------
  // build cuts  
  // ----------------------------------------------------------

  char  ch2[2048]="";
  if((cwb_lag_number==-1)&&(cwb_slag_number==-1)) {
    sprintf(ch2,"netcc[%d]>%f&&!(lag[%d]==0&&slag[%d]==0)",pp_inetcc,T_cor,nIFO,nIFO);
  } 
  if((cwb_lag_number>=0)&&(cwb_slag_number>=0)) {
    sprintf(ch2,"netcc[%d]>%f&&lag[%d]==%d&&slag[%d]==%d",
            pp_inetcc,T_cor,nIFO,cwb_lag_number,nIFO,cwb_slag_number);
  }
  if((cwb_lag_number>=0)&&(cwb_slag_number==-1)) {
    sprintf(ch2,"netcc[%d]>%f&&lag[%d]==%d",pp_inetcc,T_cor,nIFO,cwb_lag_number);
  }
  if((cwb_lag_number==-1)&&(cwb_slag_number>=0)) {
    sprintf(ch2,"netcc[%d]>%f&&slag[%d]==%d",pp_inetcc,T_cor,nIFO,cwb_slag_number);
  }

  if(simulation) sprintf(ch2,"netcc[%d]>%f",pp_inetcc,T_cor);
  else           sprintf(ch2,"%s&&(rho[%d]>%f)",ch2,pp_irho,T_cut);

  if (T_vED>0)    sprintf(ch2,"%s&&neted[0]/ecor<%f",ch2,T_vED);
  if (T_scc>0)    sprintf(ch2,"%s&&netcc[3]>%f",ch2,T_scc);
  if (T_pen>0)    sprintf(ch2,"%s&&penalty>%f",ch2,T_pen);
  if (T_hrss>0)   sprintf(ch2,"%s&&abs((TMath::Log10(hrss[%i])-TMath::Log10(hrss[%i])))<%f",ch2,i_hrss1,i_hrss2,T_hrss);

  for(int i=0;i<nFCUT;i++) {
    sprintf(ch2,"%s&&((frequency[0]<%f)||(frequency[0]>=%f))",ch2,lowFCUT[i],highFCUT[i]);
  }

/*
#if defined (CAT2_VETO) || defined (HVETO_VETO) || defined (CAT3_VETO) || defined (PEM_VETO)
  for(int i=0;i<nvdqf;i++) {
    bool belongToNet=false;
    for(int j=0;j<nIFO;j++) if(TString(vdqf[i].ifo)==ifo[j]) belongToNet=true;
    if(!belongToNet) {
      cout << "cwb_epparameters.C : Error in vdqf list defined in user_pparameters.C" << endl;
      cout << "                     it includes ifo " << vdqf[i].ifo 
           << " which is not declared in the network" << endl << endl;
      gSystem->Exit(1);  
    }
  }
#endif
*/

#ifdef CAT2_VETO
  for(int i=0;i<nvdqf;i++) {
    if((vdqf[i].cat==CWB_CAT0)||(vdqf[i].cat==CWB_CAT1)||(vdqf[i].cat==CWB_CAT2)) VDQF[nVDQF++]=vdqf[i];
  }
  char  veto_cat2[1024] = "";
  for(int i=0;i<nvdqf;i++) {
    if(vdqf[i].cat==CWB_CAT2) {
      if(strlen(veto_cat2)==0) {
        sprintf(veto_cat2,"%s_%s",CWB_CAT_NAME[vdqf[i].cat].Data(),vdqf[i].ifo);
      } else {
        sprintf(veto_cat2,"%s&&%s_%s",veto_cat2,CWB_CAT_NAME[vdqf[i].cat].Data(),vdqf[i].ifo);
      }
    }
  }
  cout << "veto_cat2 : " << veto_cat2 << endl;
  sprintf(ch2,"%s&&(%s)",ch2,veto_cat2);
#endif

#if defined (HVETO_VETO) 
  for(int i=0;i<nvdqf;i++) if(vdqf[i].cat==CWB_HVETO) VDQF[nVDQF++]=vdqf[i];

  for(int i=0;i<nvdqf;i++) {
    if(vdqf[i].cat==CWB_HVETO) {
      if(strlen(veto_not_vetoed)==0) {
        sprintf(veto_not_vetoed,"!%s_%s",CWB_CAT_NAME[vdqf[i].cat].Data(),vdqf[i].ifo);
      } else {
        char _veto_not_vetoed[1024];strcpy(_veto_not_vetoed,veto_not_vetoed);
        sprintf(veto_not_vetoed,"%s&&!%s_%s",_veto_not_vetoed,CWB_CAT_NAME[vdqf[i].cat].Data(),vdqf[i].ifo);
      }
      if(strlen(veto_vetoed)==0) {
        sprintf(veto_vetoed,"%s_%s",CWB_CAT_NAME[vdqf[i].cat].Data(),vdqf[i].ifo);
      } else {
        sprintf(veto_vetoed,"%s||%s_%s",veto_vetoed,CWB_CAT_NAME[vdqf[i].cat].Data(),vdqf[i].ifo);
      }
    }
  }
#endif

#if defined (CAT3_VETO) 
  for(int i=0;i<nvdqf;i++) if(vdqf[i].cat==CWB_CAT3) VDQF[nVDQF++]=vdqf[i];

  for(int i=0;i<nvdqf;i++) {
    if(vdqf[i].cat==CWB_CAT3) {
      if(strlen(veto_not_vetoed)==0) {
        sprintf(veto_not_vetoed,"!%s_%s",CWB_CAT_NAME[vdqf[i].cat].Data(),vdqf[i].ifo);
      } else {
        sprintf(veto_not_vetoed,"%s&&!%s_%s",veto_not_vetoed,CWB_CAT_NAME[vdqf[i].cat].Data(),vdqf[i].ifo);
      }
      if(strlen(veto_vetoed)==0) {
        sprintf(veto_vetoed,"%s_%s",CWB_CAT_NAME[vdqf[i].cat].Data(),vdqf[i].ifo);
      } else {
        sprintf(veto_vetoed,"%s||%s_%s",veto_vetoed,CWB_CAT_NAME[vdqf[i].cat].Data(),vdqf[i].ifo);
      }
    }
  }
#endif

#if defined (PEM_VETO)
  for(int i=0;i<nvdqf;i++) if(vdqf[i].cat==CWB_PEM) VDQF[nVDQF++]=vdqf[i];
#endif

#if defined (USER_VETO)
  for(int i=0;i<nvdqf;i++) if(vdqf[i].cat==CWB_USER) VDQF[nVDQF++]=vdqf[i];
#endif

  cout << "ch2 : " << ch2 << endl;

#if defined (HVETO_VETO) || defined (CAT3_VETO)
  char _veto_vetoed[1024]; strcpy(_veto_vetoed,veto_vetoed);
  sprintf(veto_vetoed,"%s&&(%s)",ch2,_veto_vetoed);
  char _veto_not_vetoed[1024]; strcpy(_veto_not_vetoed,veto_not_vetoed);
  sprintf(veto_not_vetoed,"%s&&(%s)",ch2,_veto_not_vetoed);
  cout << "veto_vetoed : " << veto_vetoed << endl;
  cout << "veto_not_vetoed : " << veto_not_vetoed << endl;
#endif
}
