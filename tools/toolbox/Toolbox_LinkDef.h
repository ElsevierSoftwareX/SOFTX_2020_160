#ifdef __CINT__ 

#pragma link off all globals; 
#pragma link off all classes; 
#pragma link off all functions; 

#pragma link C++ enum CWB_CAT+;

#pragma link C++ class CWB::Toolbox;
#pragma link C++ class CWB::CBCTool;
#pragma link C++ struct frfile+;
#pragma link C++ struct dqfile+;
#pragma link C++ struct slag+;
#pragma link C++ struct mdcshift+;
#pragma link C++ struct ifoparms+;
//#pragma link C++ struct waveSegment;
//#pragma link C++ class vector<waveSegment>+;
#pragma link C++ class vector<slag>+;
#pragma link C++ class vector<TString>+;
//#pragma link C++ class vector<double>+;
#pragma link C++ class std::map<int, TString>;
#pragma link C++ class std::pair<int,int>;
#pragma link C++ class std::map<std::pair<int,int>, int> ;

#pragma link C++ function PoissonIFunction;
#pragma link C++ function logNfit;
#pragma link C++ function DrawMDC;
#pragma link C++ function DrawWAVE;
#pragma link C++ function DrawLIVE;
#pragma link C++ function ScanMDC;
#pragma link C++ function ScanWAVE;
#pragma link C++ function ScanLIVE;
#pragma link C++ function PrintMDC;
#pragma link C++ function PrintWAVE;
#pragma link C++ function PrintLIVE;
#pragma link C++ function GetMDC;
#pragma link C++ function GetWAVE;
#pragma link C++ function GetLIVE;
#pragma link C++ function GetFileLabel;
#pragma link C++ function GetProcInfo;
#pragma link C++ function GetStart;
#pragma link C++ function GetString;
#pragma link C++ function CheckAnalysis;
#pragma link C++ function GetGitInfos;
#pragma link C++ function PrintLogoCWB;
#pragma link C++ function Draw;
#pragma link C++ function GetLiveTime;
#pragma link C++ function ReadInjType;
#pragma link C++ function MakePlotsHtmlTable;
#pragma link C++ function MakePlotsHtmlCellTable;
#pragma link C++ function GetPrecision;
#pragma link C++ function AddRho2FAR;
#pragma link C++ function GetLALVersion;

#endif // __CINT__
