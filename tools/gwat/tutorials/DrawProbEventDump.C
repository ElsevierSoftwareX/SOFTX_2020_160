//
// Read event probability skymap from output cWB text file and draw 
// Author : Gabriele Vedovato

{
  #define PROJECTION ""
  //#define PROJECTION "hammer"
  #define RESOLUTION  2
  #define COORDINATES "Geographic"

  #define EVENT_DUMP 	"eventDump.txt"

  //#define OFILE_NAME	"eventDump.png"


  gskymap* gSM = new gskymap(int(6));
  gSM->SetOptions(PROJECTION,COORDINATES,RESOLUTION/2);

  ifstream in;
  in.open(EVENT_DUMP);
  if(!in.good()) {cout << "Error Opening File : " << EVENT_DUMP << endl;gSystem->Exit(1);}

  char istring[1024];
  bool found=false;
  while(1) {
    in.getline(istring,1024);
    if (!in.good()) break;
    TObjArray* token = TString(istring).Tokenize(TString(' '));
    TObjString* stoken =(TObjString*)token->At(0);
    TString jobLabel = stoken->GetString();
    if(jobLabel.CompareTo("#skyID")==0) {found=true;continue;}

    if(found) { 
      //cout << istring << endl; 
      stoken =(TObjString*)token->At(2);
      TString sdec = stoken->GetString();
      double DEC = sdec.Atof();
      stoken =(TObjString*)token->At(5);
      TString sra = stoken->GetString();
      double RA = sra.Atof();
      //cout << DEC << " " << RA << endl;

      stoken =(TObjString*)token->At(1);
      TString sth = stoken->GetString();
      double theta = sth.Atof();
      stoken =(TObjString*)token->At(4);
      TString sph = stoken->GetString();
      double phi = sph.Atof();
      //cout << theta << " " << phi << endl;

      stoken =(TObjString*)token->At(7);
      TString sprob = stoken->GetString();
      double prob = sprob.Atof();

      int l = gSM->getSkyIndex(theta,phi);
      gSM->set(l,prob);

    }
  }

  in.close();

  TString title = TString("Probability SkyMap");

  gSM->SetTitle(title);
  gSM->Draw(0);

#ifdef OFILE_NAME
  cout << "Write : " << OFILE_NAME << endl;
  gSM->Print(OFILE_NAME);
  exit(0);
#endif

}
