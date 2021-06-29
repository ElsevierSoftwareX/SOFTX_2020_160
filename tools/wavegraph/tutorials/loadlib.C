{

  // load wavegraph library
  TString wat_path = gSystem->Getenv("HOME_WAT");
  TString wavegraph_path = wat_path+"/tools/wavegraph/lib/wavegraph.so";
  printf("Loading wavegraph          : %s ...\n",wavegraph_path.Data());
  if(gSystem->Load(wavegraph_path)) {
    cout << "error loading wavegraph library" << endl;
    gSystem->Exit(1);
  }


}
