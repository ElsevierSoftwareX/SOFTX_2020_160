{
  gskymap gsm(int(6));

  for(int i=0;i<gsm.size();i++) gsm.set(i,i);

  gsm.Draw();

  gsm.SetOptions("LVC experiment", 300,40, 1200, 670);

  gsm.Draw();
}
