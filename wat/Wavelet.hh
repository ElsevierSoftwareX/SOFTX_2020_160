/*
# Copyright (C) 2019 Sergey Klimenko
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


// Wavelet Analysis Tool
//$Id: Wavelet.hh,v 1.2 2001/12/18 18:50:13 klimenko Exp $

#ifndef WAVELET_HH
#define WAVELET_HH

#include "watfun.hh"
#include <iostream>
#include "TNamed.h"

using std::cout; 
using std::endl;

//namespace datacondAPIwat {
//namespace wat {

//: constants which rule the boundary processing
enum BORDER {B_PAD_ZERO, 
	     B_CYCLE, 
	     B_MIRROR, 
	     B_PAD_EDGE, 
	     B_POLYNOM};


//: wavelet types
enum WAVETYPE {HAAR, 
	       BIORTHOGONAL,
	       DAUBECHIES, 
	       SYMLET,
	       MEYER,
	       WDMT};

class Wavelet : public TNamed			// Wavelet base class
{
   public:

// constructors

      //: constructor
      Wavelet(int mH=1, int mL=1, int tree=0, enum BORDER border=B_CYCLE);

      //: copy constructor
      Wavelet(const Wavelet &);

      //: Virtual destructor
      virtual ~Wavelet();

      //: duplicate on heap
      //!return: Wavelet* - duplicate of *this, allocated on heap
      virtual Wavelet* Clone() const;

// access functions

      //: get array index of the first sample for (level,layer)
      virtual int getOffset(int,int);             
      //: get offset from level and frequency index
      virtual int convertF2O(int,int);
      //: get frequency index from level and offset
      virtual int convertO2F(int,int);
      //: get frequency index for (level,layer)
      virtual int convertL2F(int,int);
      //: get layer index for (level,frequency index)
      virtual int convertF2L(int,int);

// mutators

      //: set level
      inline virtual void reset() { m_Level = 0; }
      inline virtual void setLevel(int level) { m_Level = level; };
      inline virtual  int getLevel() { return m_Level; };
      inline virtual  int maxLayer() { return m_Layer; };
      inline virtual void parity(bool f) { m_Parity = f; };
      inline virtual void heterodine(bool f) { m_Heterodine = f; };
      inline virtual bool parity() { return m_Parity; };
      inline virtual bool heterodine() { return m_Heterodine; };
      inline virtual  int getMaxLevel(int);

      //: check type of wavelet tree
      inline bool BinaryTree(){ return (m_TreeType) ? true : false; }

      // print wavelet parameters
      void print();             // *MENU*
      virtual void Browse(TBrowser *b) {print();}

// data members

      //: wavelet type
      enum WAVETYPE m_WaveType;

      //: borders handling: see BORDER constants definitions above
      enum BORDER m_Border;

      //: wavelet tree type: 0-diadic, 1-binary tree
      int m_TreeType;

      //: current level of decomposition
      int m_Level;              

      //: number of layers
      int m_Layer;              

      //: number of highpass wavelet filter coefficients
      int m_H;

      //: number of lowpass wavelet filter coefficients
      int m_L;

      bool m_Heterodine;   // default is false
      bool m_Parity;       // default is true (0 delay for symmetric wavelets)

      ClassDef(Wavelet,1)

}; // class Wavelet

// inlines
inline int Wavelet::getMaxLevel(int n)
{
   int maxLevel = 0;
   for(; (n>=2*m_H) && (n>=2*m_L) && !(n&1); n/=2) maxLevel++;
   return maxLevel;
}


//}; // namespace wat
//}; // namespace datacondAPI

#endif // WAVELET_HH

















