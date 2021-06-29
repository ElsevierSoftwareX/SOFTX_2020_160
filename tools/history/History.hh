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
                          History.h  -  description 
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

#ifndef HISTORY_H
#define HISTORY_H

#include "TObject.h"
#include "TList.h"
#include "TObjString.h"
#include "TDatime.h"
#include "TTimeStamp.h"
#include <iostream>

#include "HistoryStage.hh"
#include "HistoryDefines.hh"

#define HEADER_WIDTH 80


using namespace std;

namespace CWB {

class History : public TObject  {
public:
  History(char** StageNames = NULL, int StageNumber = 0, char** TypeNames = NULL, int TypeNumber = 0, char* FilePrefix = NULL, bool HistoryModify = false);
  History(const History& History);
  ~History();

  bool StageAllowed(char* Name);
  bool StageAlreadyPresent(char* Name);
  
  bool NameAllowed(char* Name);
  bool TypeAllowed(char* Name);

  void SetStageNames(char** StageNames, int StageNumber);
  void SetTypeNames(char** TypeNames, int TypeNumber);

  void AddLog(char* Stage, char* Log, TDatime* Time = NULL);
  void AddLog(char* Stage, char* Log, int Date, int Time);

  void AddHistory(char* Stage, char* Type, char* History, TDatime* Time = NULL);
  void AddHistory(char* Stage, char* Type, char* History, int Date, int Time);

  void  SetFilePrefix(char* FilePrefix);
  char* GetFilePrefix();

  char* GetHistory(char* StageName, char* Type);
  TDatime* GetHistoryDatime(char* StageName, char* Type);
  
  int GetLogSize(char* Stage);
  char* GetLog(char* Stage, int index);

  virtual void Browse(TBrowser *b);
  void  Print();                               // *MENU*
  void  PrintSummary();                        // *MENU*
  void  DumpToTextFile(char* FileName = NULL); // *MENU*
  void  DumpToROOTFile(char* FileName = NULL); // *MENU*

  TList* GetStageNames();
  TList* GetTypeNames();

  bool SetHistoryModify(bool Replace = true);
  bool GetHistoryModify();

  char* AddStage(char* StageName);
  char* RemoveStage(char* StageName);

  char* AddType(char* TypeName);
  char* RemoveType(char* TypeName);

  char* SetStageComment(char* Stage, char* Comment = NULL);
  char* SetTypeComment(char* Stage, char* Type, char* Comment = NULL);

  char* GetStageComment(char* Stage);
  char* GetTypeComment(char* Stage, char* Type);

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
  TTimeStamp GetCreationTimeStamp(char* Stage);

  const CWB::HistoryStage* GetStage(char* Name);

  void HistoryException(int type, const char *location, const char *msgfmt, ...);

private:
  void Init();
  void Destroy();  
  bool DuplicateNames(char** NameList, int NameNumber);
  void WriteToFile(char* FileName, bool SummaryOnly = false);
  int  GetStagePosition(char* Name);
  
  TList StageNames;
  TList TypeNames;
  
  TList StageList;
  TObjString FilePrefix;

  bool HistoryModify;
  SortOrderType SortOrder;
  bool AscendingOrder;

  long CreationDate_Sec;
  long CreationDate_NSec;
  
  ClassDef(History, 3)
};

} // end namespace

#endif
