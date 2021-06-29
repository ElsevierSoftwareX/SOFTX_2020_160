{
  network* net = new network();
  for(int n=0;n<nIFO;n++) net->add(new detector(const_cast<char*>(ifo[n])));
  net->setRunID(0);		// run number is used as seed in configPlugin
 
  TString mdc_type="signal"; 
}
