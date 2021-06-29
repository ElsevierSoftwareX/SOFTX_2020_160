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
                          HistoryStage.h  -  description
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

#ifndef HISTORYSTAGE_H
#define HISTORYSTAGE_H

#include "TObject.h"
#include "TDatime.h"
#include "TList.h"
#include "TString.h"
//#include "AalException.h"
#include "HistoryDefines.hh"
#include "TTimeStamp.h"
#include "TSystem.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

namespace CWB {

class HistoryStage : public TObject  {
public: 
	HistoryStage(char** HistoryTypes = NULL, int TypeNumber = 0, char* Name = NULL, char* Comment = NULL, TDatime* Time = NULL);
  HistoryStage(const TList& HistoryTypes, char* Name = NULL, char* Comemment = NULL, TDatime* Time = NULL);
  HistoryStage(char** HistoryTypes, int TypeNumber, char* Name, char* Comment, int Date, int Time);
  HistoryStage(const TList& HistoryTypes, char* Name, char* Comment, int Date, int Time);
  HistoryStage(const HistoryStage& HistoryStage);
	~HistoryStage();

  char* SetName(char* Name);
  char* SetComment(char* Comment);
  void SetTime(TDatime* Time);
  void SetTime(int Date, int Time);
  void SetTypes(char** HistoryTypes, int TypeNumber);
  char* SetTypeComment(char* Type, char* Comment);

  char* GetName();
  char* GetComment();
  int GetDate();
  int GetTime();  
  TDatime*GetDatime();
  char* GetTypeComment(char* Type);

  void AddLog(char* LogMsg, TDatime* Time = NULL);
  void AddLog(char* LogMsg, int Date, int Time);

  void AddHistory(char* Type, char* History, char* Comment = NULL,  bool Replace = false);

  int GetHistorySize();
  int GetLogSize();
  char* GetHistoryEntry(int index);
  char* GetHistoryEntryType(int index);
  char* GetLogEntry(int index);
  int   GetLogEntryDate(int index);
  int   GetLogEntryTime(int index);
  TDatime* GetLogEntryDatime(int index);

  char* GetHistory(char* Type);

  void SortLogs(bool Ascending = true);
  
  bool TypeAllowed(char* Type);
  bool TypeAlreadyPresent(char* Type);  

  virtual void Browse(TBrowser *b);
  void  Print();
      
  bool IsSortable() const;
  int  Compare(const TObject* Obj) const;

  char* AddType(char* TypeName);
  char* RemoveType(char* TypeName);

  SortOrderType SetSortOrder(SortOrderType SortOrder);
  SortOrderType GetSortOrder();

  bool IsSortOrderInsertion();
  bool IsSortOrderDate();
  bool IsSortOrderAlphabetical();

  bool SetAscendingSortOrder();
  bool SetDescendantSortOrder();

  bool GetAscendingSortOrder();
  bool GetDescendantSortOrder();

  void Sort();

  TTimeStamp GetCreationTimeStamp();

  void HistoryStageException(int type, const char *location, const char *msgfmt, ...);
  
private:
  int Date;
  int Time;

  TList HistoryTypes;

  int NameLength;    //Name char number + 1
  char* Name;        //[NameLength]

  TList History;
  TList Logs;

  int CommentLength; //Comment char number + 1
  char* Comment;     //[CommentLength]

  SortOrderType SortOrder;
  bool AscendingOrder;

  long CreationDate_Sec;
  long CreationDate_NSec;
                        
  void Init();
  void Destroy();
  void NameSet(char* Name);
  void CommentSet(char* Comment);

  ClassDef(HistoryStage, 3)  
};

} // end namespace

#endif
