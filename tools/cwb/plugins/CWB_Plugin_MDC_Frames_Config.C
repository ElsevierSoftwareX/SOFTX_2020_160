
#include "CWB_Plugin.h"
#include <vector>

void CWB_PluginConfig() {

  //!NOISE_MDC_SIMULATION
  // Config Plugin used to setup frame files and log MDC 'on the fly' using simulation=4

  cout << "-----> CWB_Plugin_Config_Standard_InjFrames.C -> " << " gIFACTOR " << gIFACTOR << endl;

  #define nMDC 21
  #define MDC_DIR "/home/daniel.williams/data/mdc/O2"

  if(cfg->nfactor > nMDC) {
    cout << "CWB_Plugin_Config_MDC_FRAME.C - Error : number of nfactor must be < " << nMDC << endl;
    gSystem->Exit(1);
  }

  std::vector<TString>* mdcPar;
  CWB_PLUGIN_IMPORT(std::vector<TString>*,mdcPar);
  mdcPar->resize(2);

  std::vector<std::string>* channelName;
  CWB_PLUGIN_IMPORT(std::vector<std::string>*,channelName);
  channelName->resize(2);

  std::vector<std::string>* mdcType;
  CWB_PLUGIN_IMPORT(std::vector<std::string>*,mdcType);
  mdcType->resize(nMDC);


  // --------------------------------------------------------
  // define the mdcType 
  // --------------------------------------------------------
  std::vector<std::string> xmdcType(nMDC);

  (*mdcType)[0]  = "ga_0d100";                 xmdcType[0] = "ga_D0d0001";
  (*mdcType)[1]  = "ga_1d000";                 xmdcType[1] = "ga_D0d001";
  (*mdcType)[2]  = "ga_2d500";                 xmdcType[2] = "ga_D0d0025";
  (*mdcType)[3]  = "ga_4d000";                 xmdcType[3] = "ga_D0d0040";
  
  (*mdcType)[4]  = "sg_f70_q3_elliptical";     xmdcType[4] = "sg_F70Q3-elliptical";
  (*mdcType)[5]  = "sg_f235_q3_elliptical";    xmdcType[5] = "sg_F235Q3-elliptical";
  (*mdcType)[6]  = "sg_f849_q3_elliptical";    xmdcType[6] = "sg_F849Q3-elliptical";

  (*mdcType)[7]  = "sg_f70_q9_linear";         xmdcType[7] = "";
  (*mdcType)[8]  = "sg_f100_q9_linear";        xmdcType[8] = "";
  (*mdcType)[9]  = "sg_f153_q9_elliptical";    xmdcType[9] = "";
  (*mdcType)[10] = "sg_f235_q9_linear";        xmdcType[10] = "";
  (*mdcType)[11] = "sg_f361_q9_linear";        xmdcType[11] = "";
  (*mdcType)[12] = "sg_f554_q9_elliptical";    xmdcType[12] = "";
  (*mdcType)[13] = "sg_f849_q9_elliptical";    xmdcType[13] = "sg_F849Q8d9-elliptical";
  (*mdcType)[14] = "sg_f1053_q9_linear";       xmdcType[14] = "sg_F1053Q9-linear";

  (*mdcType)[15] = "sg_f70_q100_elliptical";   xmdcType[15] = "sg_F70Q100-elliptical";
  (*mdcType)[16] = "sg_f235_q100_elliptical";  xmdcType[16] = "sg_F235Q100-elliptical";
  (*mdcType)[17] = "sg_f849_q100_elliptical";  xmdcType[17] = "sg_F849Q100-elliptical";

  (*mdcType)[18] = "wnb_150d0b100d0tau0d1";    xmdcType[18] = "wnb_D0d100F100B100";
  (*mdcType)[19] = "wnb_300d0b100d0tau0d1";    xmdcType[19] = "wnb_D0d100F250B100";
  (*mdcType)[20] = "wnb_750d0b100d0tau0d1";    xmdcType[20] = "wnb_D0d100F700B100";

  // --------------------------------------------------------
  // define the mdcPar 
  // --------------------------------------------------------

  char injectionList[1024];
  char frFile[1024];

  TString type  = (*mdcType)[gIFACTOR-1].c_str();
  TString xtype = xmdcType[gIFACTOR-1].c_str();

  //sprintf(injectionList,"%s/graven/%s.xml.gz-logfile.txt",MDC_DIR,xtype.Data());
  sprintf(injectionList,"input/graven/%s.xml.gz-logfile.txt",xtype.Data());
  sprintf(frFile,"input/%s.lst",type.Data());

  (*mdcPar)[0] = injectionList;
  (*mdcPar)[1] = frFile;

  (*channelName)[0] = "L1:SCIENCE";
  (*channelName)[1] = "H1:SCIENCE";

}
