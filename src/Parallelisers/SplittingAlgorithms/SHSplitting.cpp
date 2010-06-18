/** \file SHSplitting.cpp
 *  \brief Source of spherical harmonics load splitting algorithm
 */

// System includes
//
#include <utility>
#include <queue>
#include <set>
#include <assert.h>

// External includes
//

// Class include
//
#include "Parallelisers/SplittingAlgorithms/SHSplitting.hpp"

// Project includes
//

namespace EPMDynamo {

   SHSplitting::SHSplitting(SmartSimTrunc pSimTrunc, const int nCore, const int id)
      : LoadSplitterBase(pSimTrunc, nCore, id)
   {
      assert(false);
   }

   void SHSplitting::splitRTP(int &rR0, int &rNr, ArrayI &rTh0, ArrayI &rNth, const int id) const
   {
      assert(false);
   }

   void SHSplitting::splitFDSH(ArrayI &rR0, ArrayI &rNr, ArrayI &rMs, const int id) const
   {
      assert(false);
   }

   void SHSplitting::splitSpec(ArrayI &rLs, std::vector<ArrayI> &rMs, const int id) const
   {
      assert(false);
   }

}
