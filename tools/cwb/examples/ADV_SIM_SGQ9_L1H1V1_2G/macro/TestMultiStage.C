#define USER_CONFIG           "config/user_parameters.C"

#define JOB_FILE_INIT         "data/init_931158378_44_ADV_SIM_SGQ9_L1H1V1_2G_job1.root"
#define JOB_FILE_STRAIN       "data/strain_931158378_44_ADV_SIM_SGQ9_L1H1V1_2G_job1.root"
#define JOB_FILE_CSTRAIN      "data/cstrain_931158378_44_ADV_SIM_SGQ9_L1H1V1_2G_30_job1.root"
#define JOB_FILE_COHERENCE    "data/coherence_931158378_44_ADV_SIM_SGQ9_L1H1V1_2G_30_job1.root"
#define JOB_FILE_SUPERCLUSTER "data/supercluster_931158378_44_ADV_SIM_SGQ9_L1H1V1_2G_30_job1.root"

//#define DISPLAY_ANTENNA_PATTERN

gnetwork* gNET;

void 
TestMultiStage(int jid=1, CWB_STAGE jstage=CWB_STAGE_FULL) {

  cwb2G* CWB;

  TString jfile="";

  switch(jstage) {
  case CWB_STAGE_FULL :
    break;
  case CWB_STAGE_INIT :
    break;
  case CWB_STAGE_STRAIN :
    jfile=JOB_FILE_INIT;
    break;
  case CWB_STAGE_CSTRAIN :
    jfile=JOB_FILE_STRAIN;
    break;
  case CWB_STAGE_COHERENCE :
    jfile=JOB_FILE_CSTRAIN;
    break;
  case CWB_STAGE_SUPERCLUSTER :
    jfile=JOB_FILE_COHERENCE;
    break;
  case CWB_STAGE_LIKELIHOOD :
    jfile=JOB_FILE_SUPERCLUSTER;
     break;
  default :
    cout << "STAGE  not enabled !!!" << endl;
    exit(1);
    return;
    break;
  }

  if(jstage==CWB_STAGE_FULL || jstage==CWB_STAGE_INIT) {

    // set default CWB configuration
    CWB::config CFG;
    CFG.Import("$CWB_PARAMETERS_FILE");
    CFG.Export();
    gROOT->Macro(USER_CONFIG);
    CFG.Import();
    //CFG.Print();

    CWB = new cwb2G(CFG,jstage);
    CWB::config* cfg = CWB->GetConfig();
    // set default CWB configuration
    cfg->Import(gSystem->ExpandPathName("$CWB_MACROS/cwb_inet.C"));
    CWB->SetupStage(jstage);
    //cfg->Print();
    CWB->run(jid);

  } else {

    CWB = new cwb2G(jfile,"",jstage);
    CWB->run();
  }

#ifdef DISPLAY_ANTENNA_PATTERN
  gNET = new gnetwork(*CWB->GetNetwork());
  gskymap* gSM = gNET->GetGskymap();
  gSM->SetWorldMap();
//  gNET->print(); 
  gNET->DrawAntennaPattern(3);
  gNET->DrawSitesShortLabel(kBlack);
  gNET->DrawSites(kBlack,2.0);
  gNET->DrawSitesArms(1000000,kWhite,3.0);
#else
  exit(0);
#endif

}
