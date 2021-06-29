// this macro has been used to create the default wdmXTalk catalog used for the 2G pipeline
{    
  #define OVERLAP_CATALOG "OverlapCatalog_Lev_8_16_32_64_128_256_iNu_4_Prec_10.bin"

  int layers[6] = {8, 16, 32, 64, 128, 256};
  WDM<double>* wdm[6];
  for(int i=0; i<6; i++)wdm[i] = new WDM<double>(layers[i], layers[i], 4, 10);

  monster x(wdm, 6);
  x.write(const_cast<char*>(OVERLAP_CATALOG));
} 
