{
  #define SKYMAP_FILE "data/ced_968654036_60_sg554q8d9_obj_20_job1/L1H1V1_968654042.750_968654042.750_968654042.750/probability.root"

  //
  // Draw skymap from skymap saved in a root file
  // Author : Gabriele Vedovato

  gskymap* gSM = new gskymap(TString(SKYMAP_FILE));
  gSM->SetOptions("","Geographic");
  gSM->Draw(0);
}
