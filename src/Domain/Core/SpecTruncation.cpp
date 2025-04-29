/** \file SpecTruncation.cpp
 *  \brief Source of the implementation of the spectral space truncation information
 */

// System includes
//

// External includes
//

// Class include
//
#include "Domain/Core/SpecTruncation.hpp"

// Project includes
//

namespace EPMDynamo {

   SpecTruncation::SpecTruncation(const ArrayI &ls, const std::vector<ArrayI> &ms)
      : mNl(ls.size()), mLs(ls), mMs(ms), mNm(ls.size())
   {
      // Initialise data
      this->init();
   }

   void SpecTruncation::init()
   {
      // Set number of harmonic orders per degree
      std::vector<ArrayI>::iterator it;
      int i = 0;
      for(it = this->mMs.begin(); it != this->mMs.end(); ++it, ++i)
      {
         this->mNm(i) = (*it).size();
      }
   }

}
