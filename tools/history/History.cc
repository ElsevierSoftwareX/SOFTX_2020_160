/*
# Copyright (C) 2019 Stefano Longo, Gabriele Vedovato
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/


/***************************************************************************
                          History.cc  -  description
                             -------------------
    begin                : lun set 5 2005
    copyright            : (C) 2005 by Stefano Longo
    email                : Stefano.Longo@lnl.infn.it
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "History.hh"
#include "TRandom3.h"
#include <fstream>
#include "TFile.h"
#include "TTimeStamp.h"
#include <math.h>

CWB::History::History(char** StageNames, int StageNumber, char** TypeNames, int TypeNumber, char* FilePrefix, bool HistoryModify) {
  int i;
  TTimeStamp CreationTT;  

  Init();
  CreationDate_Sec  = CreationTT.GetSec();
  CreationDate_NSec = CreationTT.GetNanoSec();
             
  if (StageNames != NULL) {
     if (DuplicateNames(StageNames, StageNumber)) HistoryException(kBreak, "CWB::History::History", "Duplicate Stage Name");                        
     for (i = 0; i < StageNumber; i++) {
       if (StageNames[i] == NULL) HistoryException(kBreak, "CWB::History::History", "StageNames[%i] is NULL", i);
       this->StageNames.AddLast(new TObjString(StageNames[i]));
     }
  }

  if (TypeNames != NULL) {
     if (DuplicateNames(TypeNames, TypeNumber)) HistoryException(kBreak, "CWB::History::History", "Duplicate Type Name");
     for (i = 0; i < TypeNumber; i++) {
       if (TypeNames[i] == NULL) HistoryException(kBreak, "CWB::History::History", "TypeNames[%i] is NULL", i);
       this->TypeNames.AddLast(new TObjString(TypeNames[i]));
     }
  }

  if (FilePrefix != NULL) this->FilePrefix.SetString(FilePrefix);
  else this->FilePrefix.SetString("");

  this->HistoryModify = HistoryModify;
}

CWB::History::History(const History& History) : TObject(History) {
  int i;

  for (i = 0; i < History.StageNames.GetSize(); i++) {
    this->StageNames.AddLast(new TObjString(*static_cast<TObjString*>(History.StageNames.At(i))));
  }
  for (i = 0; i < History.TypeNames.GetSize(); i++) {
    this->TypeNames.AddLast(new TObjString(*static_cast<TObjString*>(History.TypeNames.At(i))));
  }
  for (i = 0; i < History.StageList.GetSize(); i++) {
    this->StageList.AddLast(new HistoryStage(*static_cast<HistoryStage*>(History.StageList.At(i))));
  }

  this->FilePrefix.SetString(History.FilePrefix.GetString().Data());
  this->HistoryModify = History.HistoryModify;
  this->SortOrder = History.SortOrder;
  this->AscendingOrder = History.AscendingOrder;
  this->CreationDate_Sec = History.CreationDate_Sec;
  this->CreationDate_NSec = History.CreationDate_NSec;
}

CWB::History::~History(){
  Destroy();
}

bool
CWB::History::StageAllowed(char* Name) {
  for (int i = 0; i < StageNames.GetSize(); i++) {
    if (strcmp(Name, static_cast<TObjString*>(StageNames.At(i))->GetName()) == 0) return true;    
  }
  return false;
}

bool
CWB::History::StageAlreadyPresent(char* Name) {
  char* StrApp;
  
  for (int i = 0; i < StageList.GetSize(); i++) {
    StrApp = static_cast<HistoryStage*>(StageList.At(i))->GetName();
    if (strcmp(StrApp, Name) == 0) {
      //delete StrApp;
      return true;
    }
    //delete StrApp;
  }
  return false;
}

bool
CWB::History::NameAllowed(char* Name) {
  return StageAllowed(Name);
}

bool
CWB::History::TypeAllowed(char* Name) {
  for (int i = 0; i < TypeNames.GetSize(); i++) {
    if (strcmp(Name, static_cast<TObjString*>(TypeNames.At(i))->GetName()) == 0) return true;
  }
  return false;
}

void
CWB::History::SetStageNames(char** StageNames, int StageNumber) {
   int i;
   char* StrApp;
   HistoryStage* tmpStage;
   
   if (StageNames == NULL) HistoryException(kBreak, "CWB::History::SetStageNames", "StageNames is NULL");

   this->StageNames.Delete();

   for (i = 0; i < StageNumber; i++) {
     if (StageNames[i] == NULL) HistoryException(kBreak, "History:SetStageNames", "StageNames[%i] is NULL", i);
     this->StageNames.AddLast(new TObjString(StageNames[i]));
   }

   for (i =0; i < StageList.GetSize(); i++) {
     //StrApp = static_cast<HistoryStage*>(StageList.At(i))->GetName();
     tmpStage = static_cast<HistoryStage*>(StageList.At(i));
     StrApp = tmpStage->GetName();
     if (!NameAllowed(StrApp)) StageList.Remove(tmpStage);
     delete StrApp;
   }
}

void
CWB::History::SetTypeNames(char** TypeNames, int TypeNumber) {
  int i;
  /*
  REMOVED 03/05/2007
  char* StrApp;
  HistoryStage* tmpStage;
  */

  if (TypeNames == NULL) HistoryException(kBreak, "CWB::History::SetTypeNames", "TypeNames is NULL");
  this->TypeNames.Delete();

  for (i = 0; i < TypeNumber; i++) {
    if (TypeNames[i] == NULL) HistoryException(kBreak, "History:SetTypeNames", "TypeNames[%i] is NULL", i);
     this->TypeNames.AddLast(new TObjString(TypeNames[i]));
  }

/*  WARNING : CHANGING A TYPE ON-THE-FLY ONLY AFFECT NEW HISTORY LINES */
/*
   for (i =0; i < StageList.GetSize(); i++) {
     tmpStage = static_cast<HistoryStage*>(StageList.At(i));
     StrApp = tmpStage->GetType();
     if (!TypeAllowed(StrApp)) StageList.Remove(tmpStage);
     delete StrApp;
   }*/
  
}

void
CWB::History::AddLog(char* Stage, char* Log, TDatime* Time) {
  HistoryStage* TempStage;
  
  if (!StageAllowed(Stage)) HistoryException(kBreak, "CWB::History::AddLog", "Stage not allowed");

  if (!StageAlreadyPresent(Stage)) {
    TempStage = new HistoryStage(TypeNames, Stage, NULL, Time);
    TempStage->SetSortOrder(SortOrder);
    if (AscendingOrder) TempStage->SetAscendingSortOrder();
    else TempStage->SetDescendantSortOrder();
    StageList.AddLast(TempStage);
    //StageList.AddLast(new HistoryStage(TypeNames, Stage, NULL, Time));
  }

  HistoryStage* tmpStage = const_cast<HistoryStage*>(GetStage(Stage));
  tmpStage->AddLog(Log, Time);
}

void
CWB::History::AddLog(char* Stage, char* Log, int Date, int Time) {
  HistoryStage* TempStage;
  
  if (!StageAllowed(Stage)) HistoryException(kBreak, "CWB::History::AddLog", "Stage not allowed");

  if (!StageAlreadyPresent(Stage)) {
    TempStage = new HistoryStage(TypeNames, Stage, NULL, Date, Time);
    TempStage->SetSortOrder(SortOrder);
    if (AscendingOrder) TempStage->SetAscendingSortOrder();
    else TempStage->SetDescendantSortOrder();
    StageList.AddLast(TempStage);
    //StageList.AddLast(new HistoryStage(TypeNames, Stage, NULL, Date, Time));
  }

  HistoryStage* tmpStage = const_cast<HistoryStage*>(GetStage(Stage));
  tmpStage->AddLog(Log, Date, Time);
}

void
CWB::History::AddHistory(char* Stage, char* Type, char* History, TDatime* Time) {
  HistoryStage* TempStage;

  if (!StageAllowed(Stage)) HistoryException(kBreak, "CWB::History::AddHistory", "Stage not allowed");

  if (!StageAlreadyPresent(Stage)) {
    TempStage = new HistoryStage(TypeNames, Stage, NULL, Time);
    TempStage->SetSortOrder(SortOrder);
    if (AscendingOrder) TempStage->SetAscendingSortOrder();
    else TempStage->SetDescendantSortOrder();
    StageList.AddLast(TempStage);
    //StageList.AddLast(new HistoryStage(TypeNames, Stage, NULL, Time));
  }
                      
  HistoryStage* tmpStage = const_cast<HistoryStage*>(GetStage(Stage));
  tmpStage->AddHistory(Type, History, NULL, HistoryModify);
}

void
CWB::History::AddHistory(char* Stage, char* Type, char* History, int Date, int Time) {
  HistoryStage* TempStage;
  
  if (!StageAllowed(Stage)) HistoryException(kBreak, "CWB::History::AddHistory", "Stage not allowed");

  if (!StageAlreadyPresent(Stage)) {
    TempStage = new HistoryStage(TypeNames, Stage, NULL, Date, Time);
    TempStage->SetSortOrder(SortOrder);
    if (AscendingOrder) TempStage->SetAscendingSortOrder();
    else TempStage->SetDescendantSortOrder();
    StageList.AddLast(TempStage);
    //StageList.AddLast(new HistoryStage(TypeNames, Stage, NULL, Date, Time));
  }

  HistoryStage* tmpStage = const_cast<HistoryStage*>(GetStage(Stage));
  tmpStage->AddHistory(Type, History, NULL, HistoryModify);  
}

void
CWB::History::SetFilePrefix(char* FilePrefix) {
  this->FilePrefix.SetString(FilePrefix);
}

char*
CWB::History::GetFilePrefix() {
  return strdup(FilePrefix.GetString().Data());
}

char*
CWB::History::GetHistory(char* StageName, char* Type) {
  int i;
  char* StrApp;

  if (!NameAllowed(StageName)) HistoryException(kBreak, "CWB::History::GetHistory", "Illegal Stage Name");
  if (!TypeAllowed(Type)) HistoryException(kBreak, "CWB::History::GetHistory", "Illegal Type Name");

  for (i = 0; i < StageList.GetSize(); i++) {
    StrApp = static_cast<HistoryStage*>(StageList.At(i))->GetName();
    if (strcmp(StrApp, StageName) == 0) {
      delete StrApp;
      return static_cast<HistoryStage*>(StageList.At(i))->GetHistory(Type);
    }
    delete StrApp;
  }
  return NULL;
}

TDatime*
CWB::History::GetHistoryDatime(char* StageName, char* Type) {
  HistoryException(kBreak, "CWB::History::GetHistoryDatime", "Not implemented !!!");
  return NULL;
}

int
CWB::History::GetLogSize(char* Stage) {
  HistoryStage* tmpStage = const_cast<HistoryStage*>(GetStage(Stage));
  return tmpStage ? tmpStage->GetLogSize() : 0;
}

char*
CWB::History::GetLog(char* Stage, int index) {
  HistoryStage* tmpStage = const_cast<HistoryStage*>(GetStage(Stage));
  return tmpStage->GetLogEntry(index);
}

void
CWB::History::Browse(TBrowser *b) {
  Print();
}

void
CWB::History::Print() {
  /*
  REMOVED 03/05/2007
  int i, j, ret;
  int MaxDate, MaxTime, MaxIndex;
  char tmpStr[256];
  char* StrApp1, *StrApp2;
  TDatime* TimePtr;  

  bool *printDone;
  */

  int ret;
  char tmpStr[256];

  TRandom3 random;
  random.SetSeed(0);

  char cfg_tmp_name[256];
  sprintf(cfg_tmp_name,"/dev/shm/%f.history.cfg",fabs(random.Uniform()));
  WriteToFile(cfg_tmp_name);  

  if (getenv("IGEC_HISTORY_VIEWER") != NULL) {
    sprintf(tmpStr, "%s %s", getenv("IGEC_HISTORY_VIEWER"), cfg_tmp_name);
    ret = system(tmpStr);
    if (ret != 0) {
      sprintf(tmpStr, "vim %s", cfg_tmp_name);
      system(tmpStr);
    }
  }
  else {
    sprintf(tmpStr, "vim %s", cfg_tmp_name);
    system(tmpStr);
  }

  sprintf(tmpStr, "rm -f %s", cfg_tmp_name);
  system(tmpStr);     
}

void
CWB::History::PrintSummary() {
  /*
  REMOVED 03/05/2007
  int i, j, ret;
  int MaxDate, MaxTime, MaxIndex;
  char tmpStr[256];
  char* StrApp1, *StrApp2;
  TDatime* TimePtr;
  bool *printDone;  
  */
  
  int ret;
  char tmpStr[256];
  

  TRandom3 random;
  random.SetSeed(0);

  char cfg_tmp_name[256];
  sprintf(cfg_tmp_name,"/dev/shm/%f.history.cfg",fabs(random.Uniform()));
  WriteToFile(cfg_tmp_name, true);

  if (getenv("IGEC_HISTORY_VIEWER") != NULL) {
    sprintf(tmpStr, "%s %s", getenv("IGEC_HISTORY_VIEWER"), cfg_tmp_name);
    ret = system(tmpStr);
    if (ret != 0) {
      sprintf(tmpStr, "vim %s", cfg_tmp_name);
      system(tmpStr);
    }
  }
  else {
    sprintf(tmpStr, "vim %s", cfg_tmp_name);
    system(tmpStr);
  }

  sprintf(tmpStr, "rm -f %s", cfg_tmp_name);
  system(tmpStr);
}

void
CWB::History::DumpToTextFile(char* FileName) {
  TDatime now;
  TString file_name;

  if (FileName != NULL) file_name = TString(FileName);
  else file_name = TString(FilePrefix.GetString().Data())+TString(now.AsString())+TString(".cfg");

  cout << "Dump To File " << file_name.Data() << endl;

  WriteToFile(const_cast<char*>(file_name.Data()));
}

void
CWB::History::DumpToROOTFile(char* FileName) {
  char fname[256];
  TFile* RootFile;
  
  if (FileName != NULL) strcpy(fname, FileName);
  else {
     TDatime now;
     TString file_name = TString(FilePrefix.GetString().Data())+TString(now.AsString())+TString(".root");
     strcpy(fname, file_name.Data());     
  }
  cout << "Dump To File " << fname << endl;

  RootFile = new TFile(fname, "RECREATE");
  this->Write();
  RootFile->Close();
  delete RootFile;  
}

TList*
CWB::History::GetStageNames() {
  int i;
  TList *TmpList;

  TmpList = new TList;
  for(i = 0; i < StageNames.GetEntries(); i++) {
    TmpList->AddLast(new TObjString(StageNames.At(i)->GetName()));
  }

  return TmpList;
}

TList*
CWB::History::GetTypeNames() {
  int i;
  TList *TmpList;

  TmpList = new TList;
  for(i = 0; i < TypeNames.GetEntries(); i++) {
    TmpList->AddLast(new TObjString(TypeNames.At(i)->GetName()));
  }

  return TmpList;
}

bool
CWB::History::SetHistoryModify(bool Modify) {
  this->HistoryModify = Modify;
  return this->HistoryModify;
}

bool
CWB::History::GetHistoryModify() {
  return HistoryModify;
}

char*
CWB::History::AddStage(char* StageName) {
   if (!HistoryModify)
      HistoryException(kBreak, "CWB::History::AddStage", "History Modify not allowed");
      
   if (NameAllowed(StageName))
      HistoryException(kBreak, "CWB::History::AddStage", "Stage %s already present", StageName);

   StageNames.AddLast(new TObjString(StageName));
   
   return StageName;
}

char*
CWB::History::RemoveStage(char* StageName) {
   int i;
   HistoryStage* tmpStage;
   char* StrApp;
   TObjString* TempString;
   
   if (!HistoryModify)
      HistoryException(kBreak, "CWB::History::RemoveStage", "History Modify not allowed");

   if (!NameAllowed(StageName))
      HistoryException(kBreak, "CWB::History::RemoveStage", "Stage %s not present", StageName);

   for (i = 0; i < StageNames.GetSize(); i++) {
     TempString = static_cast<TObjString*>(StageNames.At(i));
     if (strcmp(TempString->GetName(), StageName) == 0) StageNames.Remove(TempString);
   }

   for (i = 0; i < StageList.GetSize(); i++) {
      tmpStage = static_cast<HistoryStage*>(StageList.At(i));
      StrApp = tmpStage->GetName();
      if (!NameAllowed(StrApp)) StageList.Remove(tmpStage);
      delete StrApp;
   }

   return StageName;
}

char*
CWB::History::AddType(char* TypeName) {
  if (!HistoryModify)
     HistoryException(kBreak, "Hsitory::AddType", "History modify not allowed");

  if (TypeAllowed(TypeName))
     HistoryException(kBreak, "Hsitory::AddType" , "Type %s already present", TypeName);

  TypeNames.AddLast(new TObjString(TypeName));
  
  for(int i = 0; i < StageList.GetSize(); i++) {
    static_cast<HistoryStage*>(StageList.At(i))->AddType(TypeName);
  }

  return TypeName;
}

char*
CWB::History::RemoveType(char* TypeName) {
  TObjString* TempString;
  int i;
  
  if (!HistoryModify)
     HistoryException(kBreak, "Hsitory::RemoveType", "History modify not allowed");

  if (!TypeAllowed(TypeName))
     HistoryException(kBreak, "Hsitory::RemoveType" , "Type %s not present", TypeName);

  for (i = 0; i < TypeNames.GetSize(); i++) {
    TempString = static_cast<TObjString*>(TypeNames.At(i));
    if (strcmp(TempString->GetName(), TypeName) == 0) TypeNames.Remove(TempString);
  }

  for(i = 0; i < StageList.GetSize(); i++) {
    static_cast<HistoryStage*>(StageList.At(i))->RemoveType(TypeName);
  }

  return TypeName;
}

char*
CWB::History::SetStageComment(char* Stage, char* Comment) {
  if (!StageAllowed(Stage))
     HistoryException(kBreak, "CWB::History::SetStageComment", "Stage %s not allowed", Stage);

  if (!StageAlreadyPresent(Stage))
     HistoryException(kBreak, "CWB::History::SetStageComment", "Stage %s not present yet", Stage);

  return const_cast<HistoryStage*>(GetStage(Stage))->SetComment(Comment);
}

char*
CWB::History::SetTypeComment(char* Stage, char* Type, char* Comment) {
  if (!StageAllowed(Stage))
     HistoryException(kBreak, "CWB::History::SetTypeComment", "Stage %s not allowed", Stage);

  if (!StageAlreadyPresent(Stage))
     HistoryException(kBreak, "CWB::History::SetTypeComment", "Stage %s not present yet", Stage);

  return const_cast<HistoryStage*>(GetStage(Stage))->SetTypeComment(Type, Comment);  
}

char*
CWB::History::GetStageComment(char* Stage) {
  if (!StageAllowed(Stage))
     HistoryException(kBreak, "CWB::History::GetStageComment", "Stage %s not allowed", Stage);

  if (!StageAlreadyPresent(Stage))
     HistoryException(kBreak, "CWB::History::GetStageComment", "Stage %s not present", Stage);

  return const_cast<HistoryStage*>(GetStage(Stage))->GetComment();
}

char*
CWB::History::GetTypeComment(char* Stage, char* Type) {
  if (!StageAllowed(Stage))
     HistoryException(kBreak, "CWB::History::GetTypeComment", "Stage %s not allowed", Stage);

  if (!StageAlreadyPresent(Stage))
     HistoryException(kBreak, "CWB::History::GetTypeComment", "Stage %s not present", Stage);

  return const_cast<HistoryStage*>(GetStage(Stage))->GetTypeComment(Type);
}

SortOrderType
CWB::History::SetSortOrder(SortOrderType SortOrder) {
  int i;
  
  this->SortOrder = SortOrder;
  for (i = 0; i < StageList.GetSize(); i++) {
    static_cast<HistoryStage*>(StageList.At(i))->SetSortOrder(SortOrder);
  }
  return this->SortOrder;
}

SortOrderType
CWB::History::GetSortOrder() {
  return this->SortOrder;
}

bool
CWB::History::IsSortOrderInsertion() {
  if (SortOrder == InsertionOrder) return true;
  else return false;
}

bool
CWB::History::IsSortOrderDate() {
  if (SortOrder == ElementDate) return true;
  else return false;
}

bool
CWB::History::IsSortOrderAlphabetical() {
  if (SortOrder == Alphabetical) return true;
  else return false;
}

bool
CWB::History::SetAscendingSortOrder() {
  int i;
  
  AscendingOrder = true;
  for (i = 0; i < StageList.GetSize(); i++) {
    static_cast<HistoryStage*>(StageList.At(i))->SetAscendingSortOrder();
  }
  return AscendingOrder;
}

bool
CWB::History::SetDescendantSortOrder() {
  int i;
  
  AscendingOrder = false;
  for (i = 0; i < StageList.GetSize(); i++) {
    static_cast<HistoryStage*>(StageList.At(i))->SetDescendantSortOrder();
  }
  return AscendingOrder;
}

bool
CWB::History::GetAscendingSortOrder() {
  return AscendingOrder;
}

bool
CWB::History::GetDescendantSortOrder() {
  return !AscendingOrder;
}

void
CWB::History::Sort() {
  int i;

  StageList.Sort(AscendingOrder);
  for (i = 0; i < StageList.GetSize(); i++) {
    static_cast<HistoryStage*>(StageList.At(i))->Sort();
  } 
}

TTimeStamp
CWB::History::GetCreationTimeStamp() {
  TTimeStamp CreationTT(CreationDate_Sec, CreationDate_NSec);
  
  return CreationTT;
}

TTimeStamp
CWB::History::GetCreationTimeStamp(char* Stage) {  
  if (!StageAllowed(Stage))
     HistoryException(kBreak, "CWB::History::GetCreationDate", "Stage %s not allowed", Stage);

  if (!StageAlreadyPresent(Stage))
     HistoryException(kBreak, "CWB::History::GetCreationDate", "Stage %s not present yet", Stage);

  return const_cast<HistoryStage*>(GetStage(Stage))->GetCreationTimeStamp();
}

void
CWB::History::Init() {
  HistoryModify = false;
  SortOrder = DEFAULT_SORT_ORDER;
  AscendingOrder = DEFAULT_ASCENDING;
}

void
CWB::History::Destroy() {
}

const CWB::HistoryStage*
CWB::History::GetStage(char* Name) {
  int i;
  char* StrApp;
  
  for (i = 0; i < StageList.GetSize(); i++) {
    HistoryStage* tmpStage = static_cast<HistoryStage*>(StageList.At(i));
    StrApp = tmpStage->GetName();
    if (strcmp(StrApp, Name) == 0) {
      delete StrApp;
      return tmpStage;
    }
    delete StrApp;          
  }
  return NULL;
}

bool
CWB::History::DuplicateNames(char** NameList, int NameNumber) {
   int i, j;
   
   for (i = 0; i < NameNumber; i++) {
     for (j = i + 1; j < NameNumber; j++) {
       if (strcmp(NameList[i], NameList[j]) == 0) return true;
     }
   }
   return false;
}

void
CWB::History::WriteToFile(char* FileName, bool SummaryOnly) {
  char *StrApp1, *StrApp2, tmpStr[256];
  int i, j, length;
  TDatime* TimePtr;
  
  ofstream OutFile(FileName, ios::out);

  if (OutFile.fail()) HistoryException(kBreak, "CWB::History::WriteToFile", "Error opening output file");

  Sort();
                
  for (i = 0; i < StageList.GetSize(); i++) {
    TDatime tmpTime(static_cast<HistoryStage*>(StageList.At(i))->GetDate(),
                    static_cast<HistoryStage*>(StageList.At(i))->GetTime());    

    StrApp1 = static_cast<HistoryStage*>(StageList.At(i))->GetName();
    memset(tmpStr, '*', HEADER_WIDTH);
    tmpStr[HEADER_WIDTH] = 0;
    length = strlen(StrApp1);
    tmpStr[(HEADER_WIDTH - length) / 2 - 1] = ' ';
    for (j = 0; j < length; j++) tmpStr[(HEADER_WIDTH - length) / 2 + j] = StrApp1[j];
    tmpStr[(HEADER_WIDTH + length) / 2] = ' ';
    /*
    sprintf(tmpStr, "Stage Name: %s - Date: %i/%i/%i - Time: %i:%i:%i\n", StrApp1, tmpTime.GetDay(), tmpTime.GetMonth(),
            tmpTime.GetYear(), tmpTime.GetHour(), tmpTime.GetMinute(), tmpTime.GetSecond());
    */            
    OutFile << tmpStr << endl;
    delete StrApp1;
    StrApp1 = static_cast<HistoryStage*>(StageList.At(i))->GetComment();
    //if (StrApp1 != NULL) {
    if (strlen(StrApp1) > 0) {
      OutFile << StrApp1 << endl;
      delete StrApp1;
    }
    else {
      sprintf(tmpStr, "Date: %i/%i/%i - Time: %i:%i:%i", tmpTime.GetDay(), tmpTime.GetMonth(),
            tmpTime.GetYear(), tmpTime.GetHour(), tmpTime.GetMinute(), tmpTime.GetSecond());
      OutFile << tmpStr << endl;
    }
    memset(tmpStr, '*', HEADER_WIDTH);
    tmpStr[HEADER_WIDTH + 1] = 0;
    OutFile << tmpStr << endl;
    
    if (!SummaryOnly) {
       if (static_cast<HistoryStage*>(StageList.At(i))->GetHistorySize() > 0) {
          sprintf(tmpStr, "Stage's History:");
          OutFile << tmpStr << endl;

          for (j = 0; j < static_cast<HistoryStage*>(StageList.At(i))->GetHistorySize(); j++) {
             StrApp2 = static_cast<HistoryStage*>(StageList.At(i))->GetHistoryEntryType(j);
             OutFile << StrApp2 << ": " << endl;
             delete StrApp2;
             StrApp2 = static_cast<HistoryStage*>(StageList.At(i))->GetHistoryEntry(j);
             OutFile << StrApp2 << endl;
             delete StrApp2;
          }
       }

       if (static_cast<HistoryStage*>(StageList.At(i))->GetLogSize() > 0) {
          sprintf(tmpStr, "Stage's Logs:");
          OutFile << tmpStr << endl;

          static_cast<HistoryStage*>(StageList.At(i))->SortLogs();
          for (j = 0; j < static_cast<HistoryStage*>(StageList.At(i))->GetLogSize(); j++) {
            TimePtr = static_cast<HistoryStage*>(StageList.At(i))->GetLogEntryDatime(j);
            OutFile << TimePtr->GetDay() << "/" << TimePtr->GetMonth() << "/" << TimePtr->GetYear() << " - ";
            OutFile << TimePtr->GetHour() << ":" << TimePtr->GetMinute() << ":" << TimePtr->GetSecond() << "  ";
            delete TimePtr;
            StrApp2 = static_cast<HistoryStage*>(StageList.At(i))->GetLogEntry(j);
            OutFile << StrApp2 << endl;
            delete StrApp2;
          }
       }
     }        
  }
  OutFile << endl << endl;
  OutFile.close();
}

int
CWB::History::GetStagePosition(char* Name) {
  char* StrApp;
  int   Position = -1;

  for (int i = 0; i < StageList.GetSize(); i++) {
    StrApp = static_cast<HistoryStage*>(StageList.At(i))->GetName();
    if (strcmp(StrApp, Name) == 0) {
      Position = i;
    }
  }
  return Position;
}

void CWB::History::Streamer(TBuffer &R__b)
{
   // Stream an object of class History.
   TDatime CreationDatime;
   TTimeStamp CreationTT;

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      StageNames.Streamer(R__b);
      TypeNames.Streamer(R__b);
      StageList.Streamer(R__b);
      FilePrefix.Streamer(R__b);
      if (R__v > 1) {
         R__b >> HistoryModify;
         R__b >> (Int_t&)SortOrder;
         R__b >> AscendingOrder;
         if (R__v == 1) {
            CreationDatime.Streamer(R__b);
            CreationTT.Set(CreationDatime.GetYear(), CreationDatime.GetMonth(), CreationDatime.GetDay(), CreationDatime.GetHour(), CreationDatime.GetMinute(), CreationDatime.GetSecond(), 0, true, 0);
            CreationDate_Sec  = CreationTT.GetSec();
            CreationDate_NSec = CreationTT.GetNanoSec();
         }
         else {
           R__b >> CreationDate_Sec;
           R__b >> CreationDate_NSec;
           R__b.CheckByteCount(R__s, R__c, CWB::History::IsA());
         }         
      }
      else {
        HistoryModify = false;
        SortOrder = DEFAULT_SORT_ORDER;
        AscendingOrder = DEFAULT_ASCENDING;
        CreationDate_Sec  = CreationTT.GetSec();
        CreationDate_NSec = CreationTT.GetNanoSec();
      }      
   } else {
      R__c = R__b.WriteVersion(CWB::History::IsA(), kTRUE);
      TObject::Streamer(R__b);
      StageNames.Streamer(R__b);
      TypeNames.Streamer(R__b);
      StageList.Streamer(R__b);
      FilePrefix.Streamer(R__b);
      R__b << HistoryModify;
      R__b << (Int_t)SortOrder;
      R__b << AscendingOrder;
      R__b << CreationDate_Sec;
      R__b << CreationDate_NSec;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

void
CWB::History::HistoryException(int type, const char *location, const char *msgfmt, ...) {
  cout << location << " " << msgfmt << endl;
  exit(1);
}
      
