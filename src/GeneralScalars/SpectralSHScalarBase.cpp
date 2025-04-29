/** \file SpectralSHScalarBase.cpp
 *  \brief Source of the base of spectral space data type
 */

// System includes
//

// External includes
//

// Class include
//
#include "GeneralScalars/SpectralSHScalarBase.hpp"

// Project includes
//

namespace EPMDynamo {

   SpectralSHScalarBase::SpectralSHScalarBase(const SmartTruncation pTrunc, const bool needsL0)
      : SpectralFieldBase(pTrunc), mMinL(0)
   {
      // If scalar doesn't require the L=0 mode
      if(!needsL0)
      {
         // Check if the L=0 mode is stored on this cpu
         if(pTrunc->local()->spec()->lArray()(0) == 0)
         {
            this->mMinL = 1;
         }
      }

      this->initLShells();

      // Set the not used modes to zero
      for(int l = 0; l < this->mMinL; ++l)
      {
         this->rLShell(l).setConstant(0.0);
      }
   }

   void SpectralSHScalarBase::initLShells()
   {
      int nL = this->nL();
      for(int l = 0; l < nL; ++l)
      {
         this->mLShells.push_back(SpectralSHLShell(this->nN(), this->nM(l)));
      }
   }

   void SpectralSHScalarBase::rescale(const EPMFloat scale)
   {
      // Get total number of ls
      int nL = this->nL();

      // Rescale the used modes by scale
      for(int l = this->mMinL; l < nL; ++l)
      {
         this->rLShell(l) *= scale;
      }
   }

}
