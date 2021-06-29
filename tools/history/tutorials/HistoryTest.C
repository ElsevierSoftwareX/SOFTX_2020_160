//
// Test History Class
// Author : Gabriele Vedovato


{

/* Al costruttore vanno passati nome e numero degli stages e dei tipi ammessi per ognuno dei precedenti.
   L'elenco degli stage e dei tipi attualmente deve essere un vettore di stringhe
*/

char* Stages[2];
Stages[0] = new char[256];
Stages[1] = new char[256];
char* Types[3];
Types[0] = new char[256];
Types[1] = new char[256];
Types[2] = new char[256];
strcpy(Stages[0], "FME 1");
strcpy(Stages[1], "FME 2");
strcpy(Types[0], "Tipo 1");
strcpy(Types[1], "Tipo 2");
strcpy(Types[2], "Tipo 3");

CWB::History HistoryTest(Stages, 2, Types, 3); 

/* Al solito sono presenti una schiera di metodi per verificare se uno Stage/Tipo e' ammesso
   (StageAllowed(char* Name)/TypeAllowed(char* Name)), se e' gia' presente (StageAlreadyPresent(char* Name)), etc 
   A titolo di esempio:
*/

cout << HistoryTest.StageAllowed((char*)"FME 1") << endl;
cout << HistoryTest.StageAllowed((char*)"FME 3") << endl;
cout << HistoryTest.StageAlreadyPresent((char*)"FME 1") << endl;
cout << HistoryTest.StageAlreadyPresent((char*)"FME 3") << endl;

/* L'aggiunta di informazioni di Log avviene mediante il metodo 
      void AddLog(char* Stage, char* Log, TDatime* Time = NULL)
   oppure
      void AddLog(char* Stage, char* Log, TDatime* Time = NULL)
   In cui deve essere indicato lo stage di appartenenza e il testo del Log.
   Il log viene annotato con la data specificata, altrimenti con l'ora di sistema corrispondente all'inserimento
*/


HistoryTest.AddLog((char*)"FME 1", (char*)"LOG1");
HistoryTest.AddLog((char*)"FME 1", (char*)"LOG2");
HistoryTest.AddLog((char*)"FME 2", (char*)"LOG3");

/*
   L'inserimento di voci dell'History avviene mediante il metodo
      void AddHistory(char* Stage, char* Type, char* History, TDatime* Time = NULL)
   oppure
      void AddHistory(char* Stage, char* Type, char* History, int Date, int Time)
   con le stesse modalita' del Log
*/

HistoryTest.AddHistory((char*)"FME 2", (char*)"Tipo 1", (char*)"HIST1");

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

cout << HistoryTest.GetHistory((char*)"FME 2", (char*)"Tipo 1") << endl;

}
