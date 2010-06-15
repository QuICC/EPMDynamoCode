/** \file HorizontalDomain.cpp
 *  \brief Source of the implemenation of the horizontal domain
 */

// System includes
//

// External includes
//

// Class include
//
#include "Domain/Simulation/HorizontalDomain.hpp"

// Project includes
//
#include "General/MathConstants.hpp"

namespace EPMDynamo {

   HorizontalDomain::HorizontalDomain(const int maxL, const int maxM, const int mp)
      : HorizontalTruncation(maxL, maxM, mp), mpPhGrid(new Array(this->nPh()))
   {
      // initialise the phi grid
      this->initPhiGrid();
   }

   void HorizontalDomain::initPhiGrid()
   {
      DynamoFloat dPh = 2.0*MathConstants::PI/static_cast<DynamoFloat>(this->mp());
      dPh /= static_cast<DynamoFloat>(this->nPh());

      for(int i=0; i < this->mpPhGrid->size(); ++i)
      {
         (*this->mpPhGrid)(i) = dPh*i;
      }
   }

}
