{
  //
  // Compare two cwb configurations
  // Author : Gabriele Vedovato


  CWB::config config;
  cout << config.nIFO << endl;
  config.Import("../macros/cwb_parameters.C");
  config.Import("user_parameters.C");

  config.Export();

//  config.Dump();

  CWB::config config2 = config;
  cout << config2.nIFO << endl;
  config2.bpp=0.1;
  cout << config2.nIFO << endl;

  config.Compare(config2);

  exit(0);
}
