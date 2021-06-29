void h12ascii (TH1* h, char* fname)
{
  FILE* file = fopen(fname,"w");    
  Int_t n = h->GetNbinsX();
	  
	      
  for (Int_t i=1; i<=n; i++) {
          fprintf(file, "%g %g\n",
          h->GetBinLowEdge(i)+h->GetBinWidth(i)/2,
          h->GetBinContent(i));
  }
  fclose(file);    
}
