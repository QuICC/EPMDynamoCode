/** \file BoundaryCondition.cpp
 *  \brief Source of the high lever boundary condition implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "BoundaryConditions/BoundaryCondition.hpp"

// Project includes
//

namespace EPMDynamo {

   BoundaryCondition::BoundaryCondition(const int nL, const int nX)
      : mIsHomogeneous(false), mNl(nL), mNx(nX)
   {
      this->initLHSBCValues();
   }

   void BoundaryCondition::initLHSBCValues()
   {
      for(int l = 0; l < this->nL(); ++l)
      {
         this->mLHSBCValues.push_back(Array(nX()));
      }
   }

   DynamoFloat BoundaryCondition::getLHSBC(const int i, const int l) const
   {
      return this->mLHSBCValues.at(l)(i);
   }

   const Array& BoundaryCondition::getLHSBC(const int l) const
   {
      return this->mLHSBCValues.at(l);
   }

   Array& BoundaryCondition::rLHSBCValues(const int l)
   {
      return this->mLHSBCValues.at(l);
   }
}
