//
// Test History Class
// Author : Gabriele Vedovato

{

  #define IFILE_NAME "parameters_nslag.C"


  ifstream in;
  in.open(IFILE_NAME,ios::in);
  if (!in.good()) {cout << "Error Opening File : " << IFILE_NAME << endl;exit(1);}

  char historyBuffer[10000000];

  char istringa[8192];
  int historyLength = 0;
  TString live;
  while (1) {
    in.getline(istringa,8192);
    if (!in.good()) break;
    cout << istringa << endl;
    int len = strlen(istringa);
    istringa[len]=0x0a;
    strncpy(historyBuffer+historyLength,istringa,len+1);
    historyLength += len+1;
  }
  istringa[historyLength]=0x0;

//cout << historyBuffer << endl;

/* Al costruttore vanno passati nome e numero degli stages e dei tipi ammessi per ognuno dei precedenti.
   L'elenco degli stage e dei tipi attualmente deve essere un vettore di stringhe
*/

char* Stages[3];
Stages[0] = new char[256];
Stages[1] = new char[256];
Stages[2] = new char[256];
strcpy(Stages[0], "PRODUCTION");
strcpy(Stages[1], "SIMULATION");
strcpy(Stages[2], "POST_PRODUCTION");

char* Types[2];
Types[0] = new char[256];
Types[1] = new char[256];
strcpy(Types[0], "ROOTLOGON");
strcpy(Types[1], "PARAMETERS");

CWB::History HistoryTest(Stages, 2, Types, 2); 

/* Al solito sono presenti una schiera di metodi per verificare se uno Stage/Tipo e' ammesso
   (StageAllowed(char* Name)/TypeAllowed(char* Name)), se e' gia' presente (StageAlreadyPresent(char* Name)), etc 
   A titolo di esempio:
*/
/*
cout << HistoryTest.StageAllowed("FME 1") << endl;
cout << HistoryTest.StageAllowed("FME 3") << endl;
cout << HistoryTest.StageAlreadyPresent("FME 1") << endl;
cout << HistoryTest.StageAlreadyPresent("FME 3") << endl;
*/
/* L'aggiunta di informazioni di Log avviene mediante il metodo 
      void AddLog(char* Stage, char* Log, TDatime* Time = NULL)
   oppure
      void AddLog(char* Stage, char* Log, TDatime* Time = NULL)
   In cui deve essere indicato lo stage di appartenenza e il testo del Log.
   Il log viene annotato con la data specificata, altrimenti con l'ora di sistema corrispondente all'inserimento
*/

/*
HistoryTest.AddLog("FME 1", "LOG1");
HistoryTest.AddLog("FME 1", "LOG2");
HistoryTest.AddLog("FME 2", "LOG3");
*/
/*
   L'inserimento di voci dell'History avviene mediante il metodo
      void AddHistory(char* Stage, char* Type, char* History, TDatime* Time = NULL)
   oppure
      void AddHistory(char* Stage, char* Type, char* History, int Date, int Time)
   con le stesse modalita' del Log
*/

HistoryTest.AddHistory((char*)"PRODUCTION", (char*)"PARAMETERS", historyBuffer);
HistoryTest.AddLog((char*)"PRODUCTION", (char*)"LOG1");

/*
   Per visualizzare l'History e' possibile chiamare il metodo Print (utilizzato anche Browse),
*/

HistoryTest.Print();

/*
   Altrimenti e' possibile richiamare le singole voci appartenentu ad uno stage/tipo con il metodo GetHistory.
   Richiamato con la stringa rappresentante lo stage e il tipo, questo ritornera' una copia
   della stringa corrispondente all'History oppure, se la voce non e' presente, NULL.
   Ad esempio:
*/

//cout << HistoryTest.GetHistory("FME 2", "Tipo 1") << endl;

  TFile *froot = new TFile("FileHist.root", "RECREATE");
  HistoryTest.Write();
  froot->Close();
  exit(0);
}
