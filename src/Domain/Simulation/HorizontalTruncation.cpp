/** \file HorizontalTruncation.cpp
 *  \brief Source of the implementation of the horizontal truncation information
 */

// System includes
//

// External includes
//

// Class include
//
#include "Domain/Simulation/HorizontalTruncation.hpp"

// Project includes
//
#include "Simulations/SimulationConstants.hpp"

namespace EPMDynamo {

   const EPMFloat HorizontalTruncation::msThetaDealiasing = 3.0/2.0;

   const EPMFloat HorizontalTruncation::msPhiDealiasing = 3.0;

   HorizontalTruncation::HorizontalTruncation(const int maxL, const int maxM, const int mp)
      : mMaxL(maxL),mMaxM(maxM),mMp(mp),mNl(maxL+1),mNm(maxM/mp + 1),mNh(0), mNph(static_cast<int>(nM()*msPhiDealiasing)), mNth(static_cast<int>(nL()*this->msThetaDealiasing))
   {
      // Correct grid sizes if required (for example with CSCS grid)
      this->correctGridSizes();

      // Compute the number of harmonics in total
      this->initH();
   }

   void HorizontalTruncation::correctGridSizes()
   {
      if(SimulationConstants::isCSCSGrid())
      {
         this->mNth += 2;
      }
   }

   void HorizontalTruncation::initH()
   {
      this->mNh = 0;

      for(int k =0; k < this->nL(); ++k)
      {
         mNh += this->nM(k);
      }
   }

   ArrayI HorizontalTruncation::lArray() const
   {
      ArrayI params(this->nL());

      for(int k=0; k < this->nL(); ++k)
      {
         params(k) = k;
      }

      return params;
   }

   ArrayI HorizontalTruncation::nLArray() const
   {
      ArrayI sizes(this->nM());

      for(int k=0; k < this->nM(); ++k)
      {
         sizes(k) = this->nL() - k*this->mp();
      }

      return sizes;
   }

   ArrayI HorizontalTruncation::mArray() const
   {
      ArrayI params(this->nM());

      for(int k=0; k < this->nM(); ++k)
      {
         params(k) = k*this->mp();
      }

      return params;
   }

}
