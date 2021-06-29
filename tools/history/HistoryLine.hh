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
                          HistoryLine.h  -  description
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

#ifndef AULHISTORYLINE_H
#define AULHISTORYLINE_H

#include "TObject.h"
#include "HistoryDefines.hh"
#include "TTimeStamp.h"
#include "TBuffer.h"
#include <stdlib.h>

/**
  *@author Stefano Longo
  */

using namespace std;

namespace CWB {

class HistoryLine : public TObject  {
public: 
  HistoryLine(char* Type = NULL, char* Comment = NULL, char* History = NULL);
  HistoryLine(const HistoryLine& HistoryLine);
  ~HistoryLine();

  void SetHistory(char* Type, char* Comment, char* History);
  char* SetHistoryType(char* Type);
  char* SetHistoryComment(char* Comment);
  char* SetHistoryStr(char* History);

  char* GetHistoryType();
  char* GetHistoryComment();
  char* GetHistoryStr();

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
 
  void HistoryLineException(int type, const char *location, const char *msgfmt, ...);
 
private:
   void Init();
   void Destroy();
   void TypeSet(char* Type);
   void CommentSet(char* Comment);
   void HistorySet(char* History);
   
   int TypeLength;    //Type's char number + 1
   int HistoryLength; //History's char number + 1   

   char* Type;       //[TypeLength]
   char* History;    //[HistoryLength]
   
   int CommentLength; //Comment's char number + 1
   char* Comment;    //[CommentLength]

   SortOrderType SortOrder;
   bool AscendingOrder;

   long CreationDate_Sec;
   long CreationDate_NSec;
      
   ClassDef(HistoryLine, 3)   
};

} // end namespace

#endif
