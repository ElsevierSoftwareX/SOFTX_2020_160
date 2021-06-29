{
  #define N_IFO 4

  detectorParams dP[N_IFO] = {
                              {"L1", 30.5629,  -90.7742, 0.0, 0, ( +90-197.716),    0, (    -197.716 )},  // L1

                              {"H2", 46.4551,  -119.408, 0.0, 0, ( +90-125.999),    0, ( +45-125.999)},   // H2 LCI
                              {"H3", 46.4551,  -119.408, 0.0, 0, ( +45-125.999),    0, (    -125.999)},   // H3 LCI

                              {"V1", 43.6314,   10.5045, 0.0, 0, ( +90-70.5675),    0, (    -70.5675)},   // V1
                            };


  //TString ifo[N_IFO]={"L1","H1","V1","J1"};
  TString ifo[N_IFO]={"L1","","","J1"};
  gnetwork gNET(N_IFO,ifo,dP);

  gNET.GetGskymap()->SetWorldMap();

  gNET.DrawAntennaPattern();
  gNET.DrawAntennaPattern(2);

  gNET.DrawSitesShortLabel(kBlack);
  gNET.DrawSites(kBlack,2.0);
  gNET.DrawSitesArms(1000000,kWhite,3.0);
}
