{ 
  #define OVERLAP_CATALOG "OverlapCatalog_Lev_32_64_128_256_512_1024_iNu_4_Prec_10.bin"
  #define nRES 6

  // define resolutions of interest:   
  int layers[nRES] = {32, 64, 128, 256, 512, 1024};
  
  // create corresponding WDM transforms:
  WDM<double>* wdm[nRES];
  for(int i=0; i<nRES; i++)wdm[i] = new WDM<double>(layers[i], layers[i], 4, 10);

  // create the catalog:
  monster x(wdm, nRES);

  // write the catalog in a file for future use:
  x.write(const_cast<char*>(OVERLAP_CATALOG));
} 
