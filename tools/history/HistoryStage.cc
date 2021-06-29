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
                          HistoryStage.cpp  -  description
                             -------------------
    begin                : ven set 2 2005
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

#include "HistoryStage.hh"
#include "HistoryLine.hh"
#include "HistoryLogLine.hh"
#include "TObjString.h"

CWB::HistoryStage::HistoryStage(char** HistoryTypes, int TypeNumber, char* Name, char* Comment, TDatime* Time) {
  TTimeStamp CreationTT;     
  
  Init();
  CreationDate_Sec  = CreationTT.GetSec();
  CreationDate_NSec = CreationTT.GetNanoSec();
  if ((HistoryTypes != NULL) && (TypeNumber > 0)){
     for (int i = 0; i < TypeNumber; i++) {
       if (HistoryTypes[i] == NULL) HistoryStageException(kBreak, "CWB::HistoryStage::HistoryStage", "HistoryTypes[%i] is NULL", i);       
       this->HistoryTypes.AddLast(new TObjString(HistoryTypes[i]));
     }
  }

  if (Name != NULL) SetName(Name);
  if (Comment != NULL) SetComment(Comment);
  if (Time != NULL) {
    this->Date = Time->GetDate();
    this->Time = Time->GetTime();
  }
  else {
    TDatime tmpTime;
    this->Date = tmpTime.GetDate();
    this->Time = tmpTime.GetTime();
  }
} 

CWB::HistoryStage::HistoryStage(const TList& HistoryTypes, char* Name, char* Comment, TDatime* Time) {
  TTimeStamp CreationTT;

  Init();
  CreationDate_Sec  = CreationTT.GetSec();
  CreationDate_NSec = CreationTT.GetNanoSec();
  if (Name == NULL) HistoryStageException(kBreak, "CWB::HistoryStage::HistoryStage", "Name is NULL");

  for (int i = 0; i < HistoryTypes.GetSize(); i++) {
     TObjString *tmpString = static_cast<TObjString*>(HistoryTypes.At(i));
     this->HistoryTypes.AddLast(tmpString);
  }

  if (Name != NULL) SetName(Name);
  if (Comment != NULL) SetComment(Comment);
  if (Time != NULL) {
    this->Date = Time->GetDate();
    this->Time = Time->GetTime();
  }                                        
  else {
    TDatime tmpTime;             
    this->Date = tmpTime.GetDate();
    this->Time = tmpTime.GetTime();
  }                              
}

CWB::HistoryStage::HistoryStage(char** HistoryTypes, int TypeNumber, char* Name, char* Comment, int Date, int Time) {
  TTimeStamp CreationTT;     
  
  Init();
  CreationDate_Sec  = CreationTT.GetSec();
  CreationDate_NSec = CreationTT.GetNanoSec();
  if (HistoryTypes == NULL) HistoryStageException(kBreak, "CWB::HistoryStage::HistoryStage", "HistoryType is NULL");
  if (Name == NULL) HistoryStageException(kBreak, "CWB::HistoryStage::HistoryStage", "Name is NULL");

  for (int i = 0; i < TypeNumber; i++) {
    if (HistoryTypes[i] == NULL) HistoryStageException(kBreak, "CWB::HistoryStage::HistoryStage", "HistoryTypes[%i] is NULL", i);    
    this->HistoryTypes.AddLast(new TObjString(HistoryTypes[i]));
  }
  SetName(Name);
  SetComment(Comment);
  this->Date = Date;
  this->Time = Time;
}

CWB::HistoryStage::HistoryStage(const TList& HistoryTypes, char* Name, char* Comment, int Date, int Time) {
  TTimeStamp CreationTT;   
  
  CreationDate_Sec  = CreationTT.GetSec();
  CreationDate_NSec = CreationTT.GetNanoSec();
  if (Name == NULL) HistoryStageException(kBreak, "HistoryStage::HistoryStage", "Name is NULL");

  for (int i = 0; i < HistoryTypes.GetSize(); i++) {
     TObjString *tmpString = static_cast<TObjString*>(HistoryTypes.At(i));
     this->HistoryTypes.AddLast(tmpString);
  }

  SetName(Name);
  SetComment(Comment);
  this->Date = Date;
  this->Time = Time;
}

CWB::HistoryStage::HistoryStage(const HistoryStage& HistoryStage) : TObject(HistoryStage) {
  int i;                  
                         
  this->Date = HistoryStage.Date;
  this->Time = HistoryStage.Time;
  this->NameLength = HistoryStage.NameLength;
  if (HistoryStage.Name != NULL) this->Name = strdup(HistoryStage.Name);
  else this->Name = NULL;
  this->CommentLength = HistoryStage.CommentLength;
  if (HistoryStage.Comment != NULL) this->Comment = strdup(HistoryStage.Comment);
  else this->Comment = NULL;
                                   
  for (i = 0; i < HistoryStage.HistoryTypes.GetSize(); i++) {
    this->HistoryTypes.AddLast(new TObjString(*static_cast<TObjString*>(HistoryStage.HistoryTypes.At(i))));
  }                                
  for (i = 0; i < HistoryStage.History.GetSize(); i++) {
    this->History.AddLast(new HistoryLine(*static_cast<HistoryLine*>(HistoryStage.History.At(i))));
  }                                
  for (i = 0; i < HistoryStage.Logs.GetSize(); i++) {
    this->Logs.AddLast(new HistoryLogLine(*static_cast<HistoryLogLine*>(HistoryStage.Logs.At(i))));
  }                                
  this->SortOrder = HistoryStage.SortOrder;
  this->AscendingOrder = HistoryStage.AscendingOrder;
  this->CreationDate_Sec  = HistoryStage.CreationDate_Sec;
  this->CreationDate_NSec = HistoryStage.CreationDate_NSec;
}

CWB::HistoryStage::~HistoryStage(){
  Destroy();
}

char*
CWB::HistoryStage::SetName(char* Name) {
  if (Name == NULL) HistoryStageException(kBreak, "HistoryStage::SetName", "Name is NULL");
  NameSet(Name);
  
  return Name;
}

char*
CWB::HistoryStage::SetComment(char* Comment) {
  if (Comment == NULL) HistoryStageException(kBreak, "HistoryStage::SetComment", "Comment is NULL");
  CommentSet(Comment);
  
  return Comment;
}

void
CWB::HistoryStage::SetTime(TDatime* Time) {
  this->Date = Time->GetDate();
  this->Time = Time->GetTime();
}

void
CWB::HistoryStage::SetTime(int Date, int Time) {
  this->Date = Date;
  this->Time = Time;
}

void
CWB::HistoryStage::SetTypes(char** HistoryTypes, int TypeNumber) {
   int i;
   char* StrApp;
   HistoryLine* HistoryLineApp;
   
   if (HistoryTypes == NULL) HistoryStageException(kBreak, "HistoryStage::SetTypes", "HistoryType is NULL");

   this->HistoryTypes.Delete();

   for (i = 0; i < TypeNumber; i++) {
     this->HistoryTypes.AddLast(new TObjString(HistoryTypes[i]));
   }

   for (i = 0; i < History.GetSize(); i++) {
      HistoryLineApp = static_cast<HistoryLine*>(History.At(i));
      StrApp = HistoryLineApp->GetHistoryType();
      if (!TypeAllowed(StrApp)) History.Remove(HistoryLineApp);
      delete StrApp;
   }
}

char*
CWB::HistoryStage::SetTypeComment(char* Type, char* Comment) {
  int i;
  HistoryLine* TempHistory = NULL;
  char *TempString = NULL;

  if (!TypeAllowed(Type))
     HistoryStageException(kBreak, "CWB::HistoryStage::SetTypeComment", "Type %s not allowed", Type);

  if (!TypeAlreadyPresent(Type))
     HistoryStageException(kBreak, "CWB::HistoryStage::SetTypeComment", "Type %s not present yet", Type);

  for(i = 0; i < History.GetSize(); i++) {
    TempHistory = static_cast<HistoryLine*>(History.At(i));
    TempString = TempHistory->GetHistoryType();
    if (strcmp(TempString, Type) == 0) TempHistory->SetHistoryComment(Comment);
    delete TempString;
  }

  return Comment;
}

char*
CWB::HistoryStage::GetName() {
  if (Name == NULL) return NULL;
  else return strdup(Name);
}

char*
CWB::HistoryStage::GetComment() {
  if (Comment == NULL) return NULL;
  else return strdup(Comment);
}

int
CWB::HistoryStage::GetDate() {
  return Date;
}

int
CWB::HistoryStage::GetTime() {
  return Time;
}

TDatime*
CWB::HistoryStage::GetDatime() {
  return new TDatime(Date, Time);
}

char*
CWB::HistoryStage::GetTypeComment(char* Type) {
  int i;
  HistoryLine* TempHistory = NULL;
  char *TempString = NULL, *TempComment = NULL;
  
  if (!TypeAllowed(Type))
     HistoryStageException(kBreak, "CWB::HistoryStage::GetTypeComment", "Type %s not allowed", Type);

  if (!TypeAlreadyPresent(Type))
     HistoryStageException(kBreak, "CWB::HistoryStage::GetTypeComment", "Type %s not present yet", Type);

  for(i = 0; i < History.GetSize(); i++) {
    TempHistory = static_cast<HistoryLine*>(History.At(i));
    TempString = TempHistory->GetHistoryType();
    if (strcmp(TempString, Type) == 0) TempComment = TempHistory->GetHistoryComment();
    delete TempString;
  }

  return TempComment;
}

void
CWB::HistoryStage::AddLog(char* LogMsg, TDatime* Time) {
  HistoryLogLine* TempLog;

  TempLog = new HistoryLogLine(LogMsg, Time);
  TempLog->SetSortOrder(SortOrder);
  if (AscendingOrder) TempLog->SetAscendingSortOrder();
  else TempLog->SetDescendantSortOrder();
  Logs.AddLast(TempLog);
}

void
CWB::HistoryStage::AddLog(char* LogMsg, int Date, int Time) {
  HistoryLogLine* TempLog;

  TempLog = new HistoryLogLine(LogMsg, Date, Time);
  TempLog->SetSortOrder(SortOrder);
  if (AscendingOrder) TempLog->SetAscendingSortOrder();
  else TempLog->SetDescendantSortOrder();

  Logs.AddLast(TempLog);
}

void
CWB::HistoryStage::AddHistory(char* Type, char* History, char* Comment, bool Replace) {
  HistoryLine* TempHistory;
  char* StrApp;
  int   i;
  
  if (!TypeAllowed(Type)) HistoryStageException(kBreak, "CWB::HistoryStage::AddHistory", "Type not allowed");
  if (TypeAlreadyPresent(Type)) {
    if (!Replace) HistoryStageException(kBreak, "CWB::HistoryStage::AddHistory", "Type already present");
    else {
      for (i = 0; i < this->History.GetSize(); i++) {
         StrApp = static_cast<HistoryLine*>(this->History.At(i))->GetHistoryType();
         if (strcmp(StrApp, Type) == 0) {
            this->History.Remove(this->History.At(i));
         }
         delete StrApp;
      }
    }
  }

  TempHistory = new HistoryLine(Type, Comment, History);
  TempHistory->SetSortOrder(SortOrder);
  if (AscendingOrder) TempHistory->SetAscendingSortOrder();
  else TempHistory->SetDescendantSortOrder();
  this->History.AddLast(TempHistory);
}

int
CWB::HistoryStage::GetHistorySize() {
  return History.GetSize();
}

int
CWB::HistoryStage::GetLogSize() {
  return Logs.GetSize();
}
      
char*
CWB::HistoryStage::GetHistoryEntry(int index) {
  if ((index < 0) || (index >= History.GetSize())) HistoryStageException(kBreak, "CWB::HistoryStage::GetHistoryEntry", "Index (%i) out of bounds", index);
  else return static_cast<HistoryLine*>(History.At(index))->GetHistoryStr();
  return NULL;
}       

char*
CWB::HistoryStage::GetHistoryEntryType(int index) {
  if ((index < 0) || (index >= History.GetSize())) HistoryStageException(kBreak, "CWB::HistoryStage::GetHistoryEntryType", "Index (%i) out of bounds", index);
  else return static_cast<HistoryLine*>(History.At(index))->GetHistoryType();
  return NULL;
}

char*
CWB::HistoryStage::GetLogEntry(int index) {
  if ((index < 0) || (index >= Logs.GetSize())) HistoryStageException(kBreak, "CWB::HistoryStage::GetLogEntry", "Index (%i) out of bounds", index);
  else return static_cast<HistoryLogLine*>(Logs.At(index))->GetLogStr();
  return NULL;
}

int
CWB::HistoryStage::GetLogEntryDate(int index) {
  if ((index < 0) || (index >= Logs.GetSize())) HistoryStageException(kBreak, "CWB::HistoryStage::GetLogEntryDate", "Index (%i) out of bounds", index);
  else return static_cast<HistoryLogLine*>(Logs.At(index))->GetLogDate();
  return 0;
}

int
CWB::HistoryStage::GetLogEntryTime(int index) {
  if ((index < 0) || (index >= Logs.GetSize())) HistoryStageException(kBreak, "CWB::HistoryStage::GetLogEntryTime", "Index (%i) out of bounds", index);
  else return static_cast<HistoryLogLine*>(Logs.At(index))->GetLogTime();
  return 0;
}

TDatime*
CWB::HistoryStage::GetLogEntryDatime(int index) {
  if ((index < 0) || (index >= Logs.GetSize())) HistoryStageException(kBreak, "CWB::HistoryStage::GetLogEntryDatime", "Index (%i) out of bounds", index);
  else return static_cast<HistoryLogLine*>(Logs.At(index))->GetLogDatime();
  return NULL;
}

char*
CWB::HistoryStage::GetHistory(char* Type) {
   int i;
   char* StrApp;

   if (!TypeAllowed(Type)) HistoryStageException(kBreak, "CWB::HistoryStage::GetHistory", " Illegal Type");

   for (i = 0; i < History.GetSize(); i++) {
      StrApp = static_cast<HistoryLine*>(History.At(i))->GetHistoryType();
      if (strcmp(StrApp, Type) == 0) {
        delete StrApp;
        return static_cast<HistoryLine*>(History.At(i))->GetHistoryStr();
      }
      delete StrApp;
   }
   return NULL;
}

void
CWB::HistoryStage::SortLogs(bool Ascending) {
  Logs.Sort(Ascending);
}
  
bool
CWB::HistoryStage::TypeAlreadyPresent(char* Type) {
  char* StrApp;
  
  //for (int i = 0; i < TypeNumber; i++) {
  for (int i = 0; i < History.GetSize(); i++) {
    StrApp = static_cast<HistoryLine*>(History.At(i))->GetHistoryType();
    if (strcmp(StrApp, Type) == 0) {
      //delete StrApp;
      return true;
    }
    //delete StrApp;
  }
  return false;
}

bool
CWB::HistoryStage::TypeAllowed(char* Type) {
  //for (int i = 0; i < TypeNumber; i++) {
  for (int i = 0; i < HistoryTypes.GetSize(); i++) {
    //if (strcmp(Type, HistoryTypes[i]) == 0) return true;
    if (strcmp(Type, static_cast<TObjString*>(HistoryTypes.At(i))->GetName()) == 0) return true;
  }
  return false;
}

void
CWB::HistoryStage::Browse(TBrowser *b) {
  Print();
}

void
CWB::HistoryStage::Print() {
  TDatime tmpTime(Date, Time);
  int i;
  
  cout << "Stage Time : " <<  tmpTime.GetYear() << "/" << tmpTime.GetMonth() << "/" << tmpTime.GetDay() << " - ";
  cout << tmpTime.GetHour() << ":" << tmpTime.GetMinute() << ":" << tmpTime.GetSecond() << endl;
  cout << Name << "'s History" << endl;

  for (i = 0; i < History.GetSize(); i++) {
    static_cast<HistoryLine*>(History.At(i))->Print();
  }

  cout << Name << "'s Logs" << endl;
  for (i = 0; i < Logs.GetSize(); i++) {
    static_cast<HistoryLogLine*>(Logs.At(i))->Print();
  }
}
      
bool
CWB::HistoryStage::IsSortable() const{
  return true;
}

int
CWB::HistoryStage::Compare(const TObject* Obj) const{
  char* StrTmp = NULL;
  int Result;
             
  switch(SortOrder) {
    case InsertionOrder : 
       if (this->CreationDate_Sec < static_cast<HistoryStage*>(const_cast<TObject*>(Obj))->CreationDate_Sec) Result = -1;
       else if (this->CreationDate_Sec > static_cast<HistoryStage*>(const_cast<TObject*>(Obj))->CreationDate_Sec) Result = 1;
       else if (this->CreationDate_NSec < static_cast<HistoryStage*>(const_cast<TObject*>(Obj))->CreationDate_NSec) Result = -1;
            else if (this->CreationDate_NSec > static_cast<HistoryStage*>(const_cast<TObject*>(Obj))->CreationDate_NSec) Result = 1;
                 else Result = 0;
       break;
    case ElementDate :
       if (this->Date < static_cast<HistoryStage*>(const_cast<TObject*>(Obj))->Date) Result = -1;
       else if (this->Date > static_cast<HistoryStage*>(const_cast<TObject*>(Obj))->Date) Result = 1;
       else  Result = 0;
       break;
    case Alphabetical  :
       StrTmp = static_cast<HistoryStage*>(const_cast<TObject*>(Obj))->GetName();
       Result = strcmp(this->Name, StrTmp);
       delete StrTmp;
       if (Result == 0) {
         StrTmp = static_cast<HistoryStage*>(const_cast<TObject*>(Obj))->GetComment();
         Result = strcmp(this->Comment, StrTmp);
         delete StrTmp;
       }
       break;
    default            :
//       HistoryStageException(kBreak, "CWB::HistoryStage::Compare", "Sort order not supported");
       exit(1);
       break;
  }

  return Result;
}

char*
CWB::HistoryStage::AddType(char* TypeName) {
  if (TypeAllowed(TypeName))
     HistoryStageException(kBreak, "CWB::HistoryStage::AddType", "Type %s already present", TypeName);

  HistoryTypes.AddLast(new TObjString(TypeName));

  return TypeName;
}

char*
CWB::HistoryStage::RemoveType(char* TypeName) {
   HistoryLine* tmpHistory;
   TObjString* TempString;
   char* StrApp;
   int i;
   
   if (!TypeAllowed(TypeName))
      HistoryStageException(kBreak, "CWB::HistoryStage::RemoveType", "Type %s not present", TypeName);

   for (i = 0; i < HistoryTypes.GetSize(); i++) {
     TempString = static_cast<TObjString*>(HistoryTypes.At(i));
     if (strcmp(TempString->GetName(), TypeName) == 0) HistoryTypes.Remove(TempString);
   }
   
   for (i = 0; i < History.GetSize(); i++) {
      tmpHistory = static_cast<HistoryLine*>(History.At(i));
      StrApp = tmpHistory->GetHistoryType();
      if (!TypeAllowed(StrApp)) History.Remove(tmpHistory);
      delete StrApp;
    }

    return TypeName;
}

SortOrderType
CWB::HistoryStage::SetSortOrder(SortOrderType SortOrder) {
  int i;
  
  this->SortOrder = SortOrder;
  for (i = 0; i < History.GetSize(); i++) {
    static_cast<HistoryLine*>(History.At(i))->SetSortOrder(SortOrder);
  }
  for (i = 0; i < Logs.GetSize(); i++) {
    static_cast<HistoryLogLine*>(Logs.At(i))->SetSortOrder(SortOrder);
  }
  return this->SortOrder;
}

SortOrderType
CWB::HistoryStage::GetSortOrder() {
  return SortOrder;
}

bool
CWB::HistoryStage::IsSortOrderInsertion() {
  if (SortOrder == InsertionOrder) return true;
  else return false;
}

bool
CWB::HistoryStage::IsSortOrderDate() {
  if (SortOrder == ElementDate) return true;
  else return false;
}

bool
CWB::HistoryStage::IsSortOrderAlphabetical() {
  if (SortOrder == Alphabetical) return true;
  else return false;
}

bool
CWB::HistoryStage::SetAscendingSortOrder() {
  int i;
   
  AscendingOrder = true;
  for (i = 0; i < History.GetSize(); i++) {
    static_cast<HistoryLine*>(History.At(i))->SetAscendingSortOrder();
  }
  for (i = 0; i < Logs.GetSize(); i++) {
    static_cast<HistoryLogLine*>(Logs.At(i))->SetAscendingSortOrder();
  }
  return AscendingOrder;
}

bool
CWB::HistoryStage::SetDescendantSortOrder() {
  int i;
  
  AscendingOrder = false;
  for (i = 0; i < History.GetSize(); i++) {
    static_cast<HistoryLine*>(History.At(i))->SetDescendantSortOrder();
  }
  for (i = 0; i < Logs.GetSize(); i++) {
    static_cast<HistoryLogLine*>(Logs.At(i))->SetDescendantSortOrder();
  }  
  return AscendingOrder;
}

bool
CWB::HistoryStage::GetAscendingSortOrder() {
  return AscendingOrder;
}

bool
CWB::HistoryStage::GetDescendantSortOrder() {
  return !AscendingOrder;
}

void
CWB::HistoryStage::Sort() {
   History.Sort(AscendingOrder);
   Logs.Sort(AscendingOrder);
}

TTimeStamp
CWB::HistoryStage::GetCreationTimeStamp() {
  TTimeStamp CreationTT(CreationDate_Sec, CreationDate_NSec);
  return CreationTT;
}

void
CWB::HistoryStage::Init() {
  NameLength = 0;
  Name = NULL;
  
  /*
  CommentLength = 0;
  Comment = NULL;
  */
  
  CommentLength = 1;  //Da correggere: non serializzare stringhe nulle, serializzare stringhe [0]=0
  Comment = new char[1];
  Comment[0]=0;
  SortOrder = DEFAULT_SORT_ORDER;
  AscendingOrder = DEFAULT_ASCENDING;  
}

void
CWB::HistoryStage::Destroy() {
  if (Name != NULL) delete Name;
  if (Comment != NULL) delete Name;
}

void
CWB::HistoryStage::NameSet(char* Name) {
  
   if (this->Name != NULL) delete this->Name;
  
   if (Name != NULL) {
     NameLength = strlen(Name) + 1;
     this->Name = new char[NameLength];
     strcpy(this->Name, Name);
   }
   else {
     NameLength = 1;
     this->Name = new char[1];
     this->Name[0] = 0;
   }
}

void
CWB::HistoryStage::CommentSet(char* Comment) {
  if (this->Comment != NULL) delete this->Comment;

  if (Comment != NULL) {
    CommentLength = strlen(Comment) + 1;
    this->Comment = new char[CommentLength];
    strcpy(this->Comment, Comment);
  }
  else {
    CommentLength = 1;
    this->Comment = new char[1];
    this->Comment[0] = 0;
  }
}

void CWB::HistoryStage::Streamer(TBuffer &R__b)
{
   // Stream an object of class HistoryStage.
   TDatime CreationDatime;
   TTimeStamp CreationTT;

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> Date;
      R__b >> Time;
      HistoryTypes.Streamer(R__b);
      R__b >> NameLength;
      delete [] Name;
      Name = new char[NameLength];
      R__b.ReadFastArray(Name,NameLength);
      History.Streamer(R__b);
      Logs.Streamer(R__b);
      if (R__v > 1) {
         R__b >> CommentLength;
         delete [] Comment;
         Comment = new char[CommentLength];
         R__b.ReadFastArray(Comment,CommentLength);
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
           R__b.CheckByteCount(R__s, R__c, CWB::HistoryStage::IsA());
         }         
      }
      else {
        CommentLength = 0;
        Comment = NULL;
        SortOrder = DEFAULT_SORT_ORDER;
        AscendingOrder = DEFAULT_ASCENDING;
        CreationDate_Sec  = CreationTT.GetSec();
        CreationDate_NSec = CreationTT.GetNanoSec();
      }      
   } else {
      R__c = R__b.WriteVersion(CWB::HistoryStage::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << Date;
      R__b << Time;
      HistoryTypes.Streamer(R__b);
      R__b << NameLength;
      R__b.WriteFastArray(Name,NameLength);
      History.Streamer(R__b);
      Logs.Streamer(R__b);
      R__b << CommentLength;
      R__b.WriteFastArray(Comment,CommentLength);
      R__b << (Int_t)SortOrder;
      R__b << AscendingOrder;
      R__b << CreationDate_Sec;
      R__b << CreationDate_NSec;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

void
CWB::HistoryStage::HistoryStageException(int type, const char *location, const char *msgfmt, ...) {
  cout << location << " " << msgfmt << endl;
  exit(1);  
}
