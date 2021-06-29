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
                          HistoryLogLine.cpp  -  description
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

#include "HistoryLogLine.hh"

CWB::HistoryLogLine::HistoryLogLine(char* LogStr, TDatime* Time) {
  TTimeStamp CreationTT;
  
  Init();
  CreationDate_Sec  = CreationTT.GetSec();
  CreationDate_NSec = CreationTT.GetNanoSec();
  if (Time != NULL) {
    this->Date = Time->GetDate();
    this->Time = Time->GetTime();
  }
  else {
    TDatime tmpTime;
    this->Date = tmpTime.GetDate();
    this->Time = tmpTime.GetTime();
  }

  if (LogStr != NULL) {
    LogLength = strlen(LogStr) + 1;
    Log = new char[LogLength];
    strcpy(Log, LogStr);
  }
  else {
    LogLength = 1;
    Log = new char[1];
    Log[0] = 0;
  }
}

CWB::HistoryLogLine::HistoryLogLine(char* LogStr, int Date, int Time) {
  TTimeStamp CreationTT;
  
  Init();
  CreationDate_Sec  = CreationTT.GetSec();
  CreationDate_NSec = CreationTT.GetNanoSec();
  this->Date = Date;
  this->Time = Time;
  if (LogStr != NULL) {
     LogLength = strlen(LogStr) + 1;
     Log = new char[LogLength];
     strcpy(Log, LogStr);
  }
  else {
    LogLength = 1;
    Log = new char[1];
    Log[0] = 0;
  }
}

CWB::HistoryLogLine::HistoryLogLine(const HistoryLogLine& LogLine) : TObject(LogLine) {
  this->Date = LogLine.Date;
  this->Time = LogLine.Time;
  this->LogLength = LogLine.LogLength;
  this->Log = strdup(LogLine.Log);
  this->SortOrder = LogLine.SortOrder;
  this->AscendingOrder = LogLine.AscendingOrder;
  this->CreationDate_Sec = LogLine.CreationDate_Sec;
  this->CreationDate_NSec = LogLine.CreationDate_NSec;
}

CWB::HistoryLogLine::~HistoryLogLine(){
  Destroy();
}

void
CWB::HistoryLogLine::SetLog(char* LogStr, TDatime* Time) {
  delete Log;

  this->Date = Time->GetDate();
  this->Time = Time->GetTime();
  if (LogStr != NULL) {
     LogLength = strlen(LogStr) + 1;
     Log = new char[LogLength];
     strcpy(Log, LogStr);
  }
  else {
    LogLength = 1;
    Log = new char[1];
    Log[0] = 0;
  }
}

void
CWB::HistoryLogLine::SetLog(char* LogStr, int Date, int Time) {
  delete Log;

  this->Date = Date;
  this->Time = Time;
  if (LogStr != NULL) {
     LogLength = strlen(LogStr) + 1;
     Log = new char[LogLength];
     strcpy(Log, LogStr);
  }
  else {
    LogLength = 1;
    Log = new char[1];
    Log[0] = 0;
  }
}

void
CWB::HistoryLogLine::SetLogTime(int Date, int Time) {
  this->Date = Date;
  this->Time = Time;
}

void
CWB::HistoryLogLine::SetLogTime(TDatime* Time) {
  this->Date = Time->GetDate();
  this->Time = Time->GetTime();
}

char*
CWB::HistoryLogLine::SetLogStr(char* Log) {
  delete this->Log;

  if (Log != NULL) {
     LogLength = strlen(Log) + 1;
     this->Log = new char[LogLength] + 1;
     strcpy(this->Log, Log);
  }
  else {
    LogLength = 1;
    Log = new char[1];
    Log[0] = 0;
  }

  return Log;
}

char*
CWB::HistoryLogLine::GetLogStr() {
  return strdup(Log);
}

int
CWB::HistoryLogLine::GetLogDate() {
  return Date;
}

int
CWB::HistoryLogLine::GetLogTime() {
  return Time;
}

TDatime*
CWB::HistoryLogLine::GetLogDatime() {
  return new TDatime(Date, Time);
}

void
CWB::HistoryLogLine::Browse(TBrowser *b) {
  Print();
}

void
CWB::HistoryLogLine::Print() {
  TDatime tmpTime(Date, Time);
  cout << "Log Time: " << tmpTime.GetDay() << "/" << tmpTime.GetMonth() << "/" << tmpTime.GetYear() << " - ";
  cout << tmpTime.GetHour() << ":" << tmpTime.GetMinute() << ":" << tmpTime.GetSecond() << endl;
  cout << Log << endl;
}
         
bool
CWB::HistoryLogLine::IsSortable() const {
   return true;
}

int
CWB::HistoryLogLine::Compare(const TObject* Obj) const {
  int Result;

  if (this->Date < static_cast<HistoryLogLine*>(const_cast<TObject*>(Obj))->Date) Result = -1;
       else if (this->Date > static_cast<HistoryLogLine*>(const_cast<TObject*>(Obj))->Date) Result = 1;
       else {
          if (this->Time < static_cast<HistoryLogLine*>(const_cast<TObject*>(Obj))->Time) Result = -1;
          else if (this->Time > static_cast<HistoryLogLine*>(const_cast<TObject*>(Obj))->Time) Result = 1;
          else Result = strcmp(this->Log, static_cast<HistoryLogLine*>(const_cast<TObject*>(Obj))->Log);
       }
  
  return Result;
}

SortOrderType
CWB::HistoryLogLine::SetSortOrder(SortOrderType SortOrder) {
  this->SortOrder = SortOrder;
  return this->SortOrder;
}

SortOrderType
CWB::HistoryLogLine::GetSortOrder() {
  return SortOrder;
}

bool
CWB::HistoryLogLine::IsSortOrderInsertion() {
  if (SortOrder == InsertionOrder) return true;
  else return false;
}

bool
CWB::HistoryLogLine::IsSortOrderDate() {
  if (SortOrder == ElementDate) return true;
  else return false;
}

bool
CWB::HistoryLogLine::IsSortOrderAlphabetical() {
  if (SortOrder == Alphabetical) return true;
  else return false;
}

bool
CWB::HistoryLogLine::SetAscendingSortOrder() {
  AscendingOrder = true;
  return AscendingOrder;
}

bool
CWB::HistoryLogLine::SetDescendantSortOrder() {
  AscendingOrder = false;
  return AscendingOrder;
}

bool
CWB::HistoryLogLine::GetAscendingSortOrder() {
  return AscendingOrder;
}

bool
CWB::HistoryLogLine::GetDescendantSortOrder() {
  return !AscendingOrder;
}

TTimeStamp
CWB::HistoryLogLine::GetCreationTimeStamp() {
  TTimeStamp CreationTT(CreationDate_Sec, CreationDate_NSec);
  return CreationTT;
}

void
CWB::HistoryLogLine::Init() {
  LogLength = 0;
  Log = NULL;
  SortOrder = DEFAULT_SORT_ORDER;
  AscendingOrder = DEFAULT_ASCENDING;
}

void
CWB::HistoryLogLine::Destroy() {
  if (Log != NULL) delete Log;
}

void CWB::HistoryLogLine::Streamer(TBuffer &R__b)
{
   // Stream an object of class HistoryLogLine.
   TDatime CreationDatime;
   TTimeStamp CreationTT;

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> Date;
      R__b >> Time;
      R__b >> LogLength;
      delete [] Log;
      Log = new char[LogLength];
      R__b.ReadFastArray(Log,LogLength);
      if (R__v > 1) {
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
         }
         R__b.CheckByteCount(R__s, R__c, CWB::HistoryLogLine::IsA());
      }
      else {
        SortOrder = DEFAULT_SORT_ORDER;
        AscendingOrder = DEFAULT_ASCENDING;
        CreationDate_Sec  = CreationTT.GetSec();
        CreationDate_NSec = CreationTT.GetNanoSec();
      }      
   } else {
      R__c = R__b.WriteVersion(CWB::HistoryLogLine::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << Date;
      R__b << Time;
      R__b << LogLength;
      R__b.WriteFastArray(Log,LogLength);
      R__b << (Int_t)SortOrder;
      R__b << AscendingOrder;
      R__b << CreationDate_Sec;
      R__b << CreationDate_NSec;
      R__b.SetByteCount(R__c, kTRUE);
   }
}
