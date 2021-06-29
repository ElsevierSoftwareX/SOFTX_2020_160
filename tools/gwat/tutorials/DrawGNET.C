{
  #define WORLD_MAP_DIR "$CWB_GWAT/data/"

  detector L1((char*)"L1");
  detector H1((char*)"H1");
  detector V1((char*)"V1");

  gnetwork gNET;
  gNET.add(&L1);
  gNET.add(&H1);
  gNET.add(&V1);

  gskymap* gSM = gNET.GetGskymap();
//  gSM->SetOptions("hammer","Geographic",2);

  TString world_map = gSystem->ExpandPathName(WORLD_MAP_DIR);
  gSM->SetWorldMapPath(world_map.Data());
  gSM->SetWorldMap();

  gNET.DrawAntennaPattern(3);

  gNET.DrawSitesShortLabel(kBlack);
  gNET.DrawSites(kBlack,2.0);
  gNET.DrawSitesArms(1000000,kWhite,3.0);

//  gNET.setSkyMaps(int(7));
//  gNET.Draw("skyMask");
}
