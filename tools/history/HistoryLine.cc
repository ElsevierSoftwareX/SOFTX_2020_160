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
                          HistoryLine.cpp  -  description
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

#include "HistoryLine.hh"
#include <iostream>
#include <string.h>
#include "TDatime.h"  //  [0.94.01]

CWB::HistoryLine::HistoryLine(char* Type, char* Comment, char* History){
   TTimeStamp CreationTT;
   
   Init();
   CreationDate_Sec  = CreationTT.GetSec();
   CreationDate_NSec = CreationTT.GetNanoSec();
   TypeSet(Type);
   CommentSet(Comment);
   HistorySet(History);
}

CWB::HistoryLine::HistoryLine(const HistoryLine& HistoryLine) : TObject(HistoryLine) {
  this->TypeLength = HistoryLine.TypeLength;
  this->CommentLength = HistoryLine.CommentLength;
  this->HistoryLength = HistoryLine.HistoryLength;
  this->Type = strdup(HistoryLine.Type);
  this->Comment = strdup(HistoryLine.Comment);
  this->History = strdup(HistoryLine.History);
  this->SortOrder = HistoryLine.SortOrder;
  this->AscendingOrder = HistoryLine.AscendingOrder;
  this->CreationDate_Sec = HistoryLine.CreationDate_Sec;
  this->CreationDate_NSec = HistoryLine.CreationDate_NSec;
}

CWB::HistoryLine::~HistoryLine(){
  Destroy();
}

void
CWB::HistoryLine::SetHistory(char* Type, char* Comment, char* History) {
  TypeSet(Type);
  CommentSet(Comment);
  HistorySet(History);
}

char*
CWB::HistoryLine::SetHistoryType(char* Type) {
  TypeSet(Type);  
  return Type;
}

char*
CWB::HistoryLine::SetHistoryComment(char* Comment) {
  CommentSet(Comment);
  return History;
}

char*
CWB::HistoryLine::SetHistoryStr(char* History) {
  HistorySet(History);
  return History;
}

char*
CWB::HistoryLine::GetHistoryType() {
  if (Type == NULL) return NULL;
  else return strdup(Type);
}

char*
CWB::HistoryLine::GetHistoryComment() {
  if (Comment == NULL) return NULL;
  else return strdup(Comment);
}

char*
CWB::HistoryLine::GetHistoryStr() {
  if (History == NULL) return NULL;
  else return strdup(History);
}

void
CWB::HistoryLine::Browse(TBrowser *b) {
  Print();
}

void
CWB::HistoryLine::Print() {
  cout << "Type    : " << Type << endl;
  cout << "Comment : " << Comment << endl;
  cout << "History : " << History << endl;
}
      
bool
CWB::HistoryLine::IsSortable() const {
   return true;
}

int
CWB::HistoryLine::Compare(const TObject* Obj) const {
  int Result;

  switch(SortOrder) {
    case InsertionOrder :
       if (this->CreationDate_Sec < static_cast<HistoryLine*>(const_cast<TObject*>(Obj))->CreationDate_Sec) Result = -1;
       else if (this->CreationDate_Sec > static_cast<HistoryLine*>(const_cast<TObject*>(Obj))->CreationDate_Sec) Result = 1;
       else if (this->CreationDate_NSec < static_cast<HistoryLine*>(const_cast<TObject*>(Obj))->CreationDate_NSec) Result = -1;
            else if (this->CreationDate_NSec > static_cast<HistoryLine*>(const_cast<TObject*>(Obj))->CreationDate_NSec) Result = 1;
            else Result = 0;
    case Alphabetical :
       Result = strcmp(this->Type, static_cast<HistoryLine*>(const_cast<TObject*>(Obj))->Type);
       if (Result == 0) Result = strcmp(this->History, static_cast<HistoryLine*>(const_cast<TObject*>(Obj))->History);
       break;
    default :
//       HistoryLineException(kBreak, "CWB::HistoryLine::Compare", "Sort Order not supported XXX");
       exit(1);
       break;
  }
  return Result;
}

SortOrderType
CWB::HistoryLine::SetSortOrder(SortOrderType SortOrder) {  
  if (SortOrder == ElementDate) this->SortOrder = InsertionOrder;   
  else this->SortOrder = SortOrder;
  return this->SortOrder;
}

SortOrderType
CWB::HistoryLine::GetSortOrder() {
  return SortOrder;
}

bool
CWB::HistoryLine::IsSortOrderInsertion() {
  if (SortOrder == InsertionOrder) return true;
  else return false;
}

bool
CWB::HistoryLine::IsSortOrderDate() {
  if (SortOrder == ElementDate) return true;
  else return false;
}

bool
CWB::HistoryLine::IsSortOrderAlphabetical() {
  if (SortOrder == Alphabetical) return true;
  else return false;
}

bool
CWB::HistoryLine::SetAscendingSortOrder() {
  AscendingOrder = true;
  return AscendingOrder;  
}

bool
CWB::HistoryLine::SetDescendantSortOrder() {
  AscendingOrder = false;
  return AscendingOrder;
}

bool
CWB::HistoryLine::GetAscendingSortOrder() {
  return AscendingOrder;
}

bool
CWB::HistoryLine::GetDescendantSortOrder() {
  return !AscendingOrder;
}

TTimeStamp
CWB::HistoryLine::GetCreationTimeStamp() {
  TTimeStamp CreationTT(CreationDate_Sec, CreationDate_NSec);
  return CreationTT;
}

void
CWB::HistoryLine::Init() {
  TypeLength = 0;
  Type = NULL;
  HistoryLength = 0;
  History = NULL;
  CommentLength = 0;
  Comment = NULL;
  SortOrder = DEFAULT_SORT_ORDER;
  AscendingOrder = DEFAULT_ASCENDING;  
}

void
CWB::HistoryLine::Destroy() {
  if (Type != NULL) delete Type;
  if (History != NULL) delete History;
  if (Comment != NULL) delete Comment;
}

void
CWB::HistoryLine::TypeSet(char* Type) {
  if (this->Type != NULL) delete this->Type;

  if (Type != NULL) {
    TypeLength = strlen(Type) + 1;
    this->Type = new char[TypeLength];
    strcpy(this->Type, Type);
  }
  else {
    TypeLength = 1;
    this->Type = new char[1];
    this->Type[0] = 0;
  }
}

void
CWB::HistoryLine::CommentSet(char* Comment) {
  if (this->Comment != NULL) delete this->Comment;

  if (Comment != NULL) {
    CommentLength = strlen(Comment) + 1;
    this->Comment = new char[CommentLength];
    strcpy(this->Comment, Comment);
  }
  else {
    HistoryLength = 1;
    this->Comment = new char[1];
    this->Comment[0] = 0;
  }
}

void
CWB::HistoryLine::HistorySet(char* History) {
   if (this->History != NULL) delete this->History;

   if (History != NULL) {
     HistoryLength = strlen(History) + 1;
     this->History = new char[HistoryLength];
     strcpy(this->History, History);
   }
   else {
     HistoryLength = 1;
     this->History = new char[1];
     this->History[0] = 0;
   }     
}

void CWB::HistoryLine::Streamer(TBuffer &R__b)
{
   // Stream an object of class HistoryLine.
   TDatime CreationDatime;
   TTimeStamp CreationTT;

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> TypeLength;
      R__b >> HistoryLength;
      delete [] Type;
      Type = new char[TypeLength];
      R__b.ReadFastArray(Type,TypeLength);
      delete [] History;
      History = new char[HistoryLength];
      R__b.ReadFastArray(History,HistoryLength);
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
            R__b.CheckByteCount(R__s, R__c, CWB::HistoryLine::IsA());
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
      R__c = R__b.WriteVersion(CWB::HistoryLine::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << TypeLength;
      R__b << HistoryLength;
      R__b.WriteFastArray(Type,TypeLength);
      R__b.WriteFastArray(History,HistoryLength);
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
CWB::HistoryLine::HistoryLineException(int type, const char *location, const char *msgfmt, ...) {
  cout << location << " " << msgfmt << endl;
  exit(1);
}

