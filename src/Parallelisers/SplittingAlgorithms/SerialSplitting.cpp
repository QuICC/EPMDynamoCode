/** \file SerialSplitting.cpp
 *  \brief Source of the implementation of a serial code "load splitting"
 */

// System includes
//

// External includes
//

// Class include
//
#include "Parallelisers/SplittingAlgorithms/SerialSplitting.hpp"

// Project includes
//

namespace EPMDynamo {

   SerialSplitting::SerialSplitting(SmartSimTrunc pSimTrunc, const int nCore, const int id)
      : LoadSplitterBase(pSimTrunc, nCore, id)
   {
   }

   void SerialSplitting::splitRTP(int &rR0, int &rNr, ArrayI &rTh0, ArrayI &rNth, const int id) const
   {
      // Set the radial truncation information
      rR0 = 0;
      rNr = sim()->rad()->nR();

      // Set the theta truncation information
      rTh0.resize(rNr);
      rTh0.setConstant(0);

      int n = sim()->hoz()->nTh();
      rNth.resize(rNr);
      rNth.setConstant(n);
   }

   void SerialSplitting::splitFDSH(ArrayI &rR0, ArrayI &rNr, ArrayI &rMs, const int id) const
   {
      // Set the radial truncation information
      rR0.resize(sim()->hoz()->nM());
      rR0.setConstant(0);

      rNr.resize(sim()->hoz()->nM());
      rNr.setConstant(sim()->rad()->nR());

      // Set the harmonic orders truncation information
      rMs = sim()->hoz()->mArray();
   }

   void SerialSplitting::splitSpec(ArrayI &rLs, std::vector<ArrayI> &rMs, const int id) const
   {
      // Set the harmonic degrees truncation information
      rLs = sim()->hoz()->lArray();

      // Set the harmonic orders truncation information
      rMs.clear();
      for(int l=0; l < sim()->hoz()->nL(); ++l)
      {
         rMs.push_back(ArrayI(sim()->hoz()->nM(l)));
         for(int j = 0; j < rMs.at(l).size(); ++j)
         {
            rMs.at(l)(j) = j;
         }
      }
   }

}
