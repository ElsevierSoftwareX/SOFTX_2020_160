/*
# Copyright (C) 2019 Gabriele Vedovato
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/


// command used to compile macros with ROOT ACLiC compiler
{

  CWB::Toolbox TB;

  TB.checkFile(gSystem->Getenv("CWB_ROOTLOGON_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_PARAMETERS_FILE"));
  TB.checkFile(gSystem->Getenv("CWB_NETC_FILE"));


  // get macro path
  TString macroPath="";
  if(gSystem->Getenv("CWB_COMPILE_MACRO_NAME")!=NULL) {
    // user define condor tag
    macroPath=TString(gSystem->Getenv("CWB_COMPILE_MACRO_NAME"));
  }
  TB.checkFile(macroPath);

  if(!macroPath.EndsWith(".C"))
    {cout << "cwb_compile - Error : macro : " << macroPath << " must have extension .C" << endl;gSystem->Exit(1);}

  TString libPath = macroPath;
  libPath.ReplaceAll(".C","_C");

  int  error;
  char command[1024];

  // get macro + cWB compiler flags
  TString macroTmp = CWB::Toolbox::addCWBFlags(TMacro(macroPath));

  // compile & load macro  
  cout << endl << "cwb_compile : Macro compilation ..." << endl << endl;
  int success = gSystem->CompileMacro(TString(macroTmp).Data(),"kf",libPath);
  if(!success) {
    cout << endl << "cwb_compile : Macro Compilation Failed !!! " << endl << endl;
    exit(1);
  }

  // check if the compiled shared files exists and is readable  
  int check = gROOT->LoadMacro((libPath+".so").Data(),&error,true);
  if(check==-1) {
    cout << endl << "cwb_compile : Macro Compilation Failed !!! " << endl << endl;
    exit(1);
  }

  // remove temporary ACLiC_dict_rdict.pcm file created by ACLiC
  TString macroTmp1 = macroPath;
  macroTmp1.ReplaceAll(".C","_C_ACLiC_dict_rdict.pcm");
  sprintf(command, "/bin/rm -f %s", macroTmp1.Data());
  error=gSystem->Exec(command);
  if(error) {
    cout << endl << "cwb_compile - Error -> " << command << endl << endl;
    exit(1);
  }
  // remove temporary _C.d file created by ACLiC
  TString macroTmp2 = macroTmp;
  macroTmp2.ReplaceAll(".C","_C.d");
  sprintf(command, "/bin/rm -f %s", macroTmp2.Data());
  error=gSystem->Exec(command);
  if(error) {
    cout << endl << "cwb_compile - Error -> " << command << endl << endl;
    exit(1);
  }

  cout << endl << "cwb_compile : Macro successfully compiled !!! " << endl << endl;
  exit(0);
}

