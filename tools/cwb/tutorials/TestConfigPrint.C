{
  //
  // Read and print cwb configuration 
  // Author : Gabriele Vedovato

  int nDQF = 0;
  dqfile DQF[20];

  CWB::config config;
  cout << config.nIFO << endl;
  config.Import("cwb_parameters.C");

  config.Export();

  nDQF=12;
  dqfile dqf[12]={
                    {"L1" ,"input/S6A_OFFLINE_L1SCIENCE.txt",         CWB_CAT0, 0., false, false},
                    {"L1" ,"input/S6A_OFFLINE_L1_DQCAT1SEGMENTS.txt", CWB_CAT1, 0., true, false},
                    {"L1" ,"input/S6A_OFFLINE_L1_DQCAT2SEGMENTS.txt", CWB_CAT2, 0., true, false},
                    {"L1" ,"input/S6A_OFFLINE_L1_DQCAT4SEGMENTS.txt", CWB_CAT1, 0., true, false},
                    {"H1" ,"input/S6A_OFFLINE_H1SCIENCE.txt",         CWB_CAT0, 0., false, false},
                    {"H1" ,"input/S6A_OFFLINE_H1_DQCAT1SEGMENTS.txt", CWB_CAT1, 0., true, false},
                    {"H1" ,"input/S6A_OFFLINE_H1_DQCAT2SEGMENTS.txt", CWB_CAT2, 0., true, false},
                    {"H1" ,"input/S6A_OFFLINE_H1_DQCAT4SEGMENTS.txt", CWB_CAT1, 0., true, false},
                    {"V1" ,"input/S6A_OFFLINE_V1SCIENCE.txt",         CWB_CAT0, 0., false, false},
                    {"V1" ,"input/S6A_OFFLINE_V1_DQCAT1SEGMENTS.txt", CWB_CAT1, 0., true, false},
                    {"V1" ,"input/S6A_OFFLINE_V1_DQCAT2SEGMENTS.txt", CWB_CAT2, 0., true, false},
                    {"V1" ,"input/S6A_OFFLINE_V1_DQCAT4SEGMENTS.txt", CWB_CAT1, 0., true, false} };
  for(int i=0;i<nDQF;i++) DQF[i]=dqf[i];

  config.Import();

//  config.Print();

  CWB::config config2 = config;
  cout << config2.nIFO << endl;
  //config2.nIFO=9;
  config2.Acore=0.1;
  config2.gamma=0.1;

  config.Compare(config2);

  config.Print();
  config.Print("test_config1.C");
/*
  config.nIFO=9;
  config.Print("test_config2.C");
*/
  exit(0);
}
