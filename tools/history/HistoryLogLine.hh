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
                          HistoryLogline.h  -  description
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

#ifndef HISTORYLOGLINE_H
#define HISTORYLOGLINE_H

#include "TObject.h"
#include "TDatime.h"
#include "HistoryDefines.hh"
#include "TTimeStamp.h"
#include "TBuffer.h"
#include <iostream>
#include <string.h>

/**
  *@author Stefano Longo
  */

using namespace std;
 
namespace CWB {

class HistoryLogLine : public TObject  {
public:
  HistoryLogLine(char* LogStr = NULL, TDatime* Time = NULL);
  HistoryLogLine(char* LogStr, int Date, int Time);
  HistoryLogLine(const HistoryLogLine& LogLine);
  ~HistoryLogLine();

  void SetLog(char* LogStr, TDatime* Time = NULL);
  void SetLog(char* LogStr, int Date, int Time);
  void SetLogTime(int Date, int Time);
  void SetLogTime(TDatime* Time);
  char* SetLogStr(char* Log);

  char* GetLogStr();
  int GetLogDate();
  int GetLogTime();
  TDatime* GetLogDatime();

  virtual void Browse(TBrowser *b);
  void  Print();
            
  bool IsSortable() const;
  int  Compare(const TObject* Obj) const;

  SortOrderType SetSortOrder(SortOrderType SortOrder);
  SortOrderType GetSortOrder();

  bool IsSortOrderInsertion();
  bool IsSortOrderDate();
  bool IsSortOrderAlphabetical();

  bool SetAscendingSortOrder();
  bool SetDescendantSortOrder();

  bool GetAscendingSortOrder();
  bool GetDescendantSortOrder();

  TTimeStamp GetCreationTimeStamp();

private:
   void Init();
   void Destroy();
   
   int Date, Time; //Date and Time in TDatetime.GetDate/Time() format;

   int LogLength;  //Lengh of log = char number + 1
   char* Log;      //[LogLength]

   SortOrderType SortOrder;
   bool AscendingOrder;

   long CreationDate_Sec;
   long CreationDate_NSec;
      
   ClassDef(HistoryLogLine, 3)   
};

} // end namespace

#endif
