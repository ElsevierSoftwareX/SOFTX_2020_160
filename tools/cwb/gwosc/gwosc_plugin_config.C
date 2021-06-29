#include "CWB_Plugin.h"

void CWB_PluginConfig() {

  // Config Plugin used to generate injections from posterior samples XML files
  // used by the cwb_gwosc command

  CWB::mdc* MDC;
  CWB_PLUGIN_IMPORT(CWB::mdc*,MDC);

  CWB::config** cfg;
  CWB_PLUGIN_IMPORT(CWB::config**,cfg);

  int* gIFACTOR;
  CWB_PLUGIN_IMPORT(int*,gIFACTOR);

  char symlink[1024];
  sprintf(symlink,"config/injections.xml");
  TString xmlFile = CWB::Toolbox::getFileName(symlink);	// get path from symbolic link
  if(xmlFile!="") {					// it is a symbolic link
    if(xmlFile[0]=='.') {				// if path is relative, start with '../' then '../' is removed 
      xmlFile = xmlFile(xmlFile.Index('/')+1, xmlFile.Sizeof()-xmlFile.Index('/')-2);
    }
    cout << endl;
    cout << "---------> CWB_PluginConfig: injections.xml = " << xmlFile << endl;
  } else {
    xmlFile=symlink;					// not a symbolic link
  }

  TString inspOptions="";
  inspOptions+= "--dir "+TString((*cfg)->nodedir)+" ";
  inspOptions+= "--xml "+xmlFile;
  MDC->SetInspiral("PE",inspOptions);
//  MDC->SetInspiralCLB(clbFile);
}

