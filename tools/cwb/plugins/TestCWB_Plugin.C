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


{
// Test
  //cwb CWB("job_data/job_968654036_sg554q8d9_obj_1_999741748.root");
  cwb CWB("config/user_parameters.C");
  CWB::config* cfg = CWB.GetConfig();
  cout << cfg->nodedir << endl;
  cfg->Import("../../../cwb/macros/cwb_inet.C");
  cout << cfg->nodedir << endl;
  int runID=1;
  //cfg->Print();

  int err=0;

  //err=CWB.SetPlugin("plugins/CWB_Plugin_Test.C","plugins/CWB_configPlugin_Test.C");

  //err=CWB.SetPlugin("plugins/CWB_Plugin_PhaseMisCal.C");

  //err=CWB.SetPlugin("plugins/CWB_Plugin_TShiftMisCal.C");

  //err=CWB.SetPlugin("plugins/CWB_Plugin_AmplitudeMisCal.C");

  //err=CWB.SetPlugin("plugins/CWB_Plugin_SimNoise.C");

  //err=CWB.SetPlugin("plugins/CWB_Plugin_MakeSpectrum.C");

  //err=CWB.SetPlugin("plugins/CWB_Plugin_InjectMDC.C");

  //cfg.cedDump=false;
  //cfg.factors[0]=0.64;
  //err=CWB.SetPlugin("plugins/CWB_Plugin_SimMDC_SimData.C");

  cfg.cedDump=true;
  cfg.factors[0]=0.64;
  //err=CWB.SetPlugin("plugins/CWB_Plugin_TestClassMDC.C");
  err=CWB.SetPlugin("plugins/CWB_Plugin_TestClassMDC.C","plugins/CWB_Plugin_TestClassMDC_Config.C");

  if(err) {cout << "Plugin Error !!!" << endl;gSystem->Exit(1);}

  CWB.run(runID);
}
