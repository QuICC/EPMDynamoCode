/** \file AssocLegendreTransformOperatorTest.cpp
 *  \brief Definitions of the AssocLegendreTransformOperatorTest methods.
 */

// System includes
//

// External includes
//

// Class include
//
#include "AssocLegendreTransformOperatorTest.hpp"
#include "Polynomials/AssociatedLegendreOperator.hpp"

// Project includes
//

namespace EPMDynamo {

   AssocLegendreTransformOperatorTest::AssocLegendreTransformOperatorTest(SmartTruncation pTrunc)
      : TransformOperatorTest<AssocLegendreBasis>(pTrunc)
   {
      // Set prefix
      this->mPrefix = "assocLeg_";

      // Set ordering string
      this->mOrderString = "_M";

      // Set max order index
      this->mOrderingN = pTrunc->sim()->hoz()->nM();

      // Set grid size
      this->mGridSize = pTrunc->sim()->hoz()->nTh();

      // Set spectrum size
      this->mSpectrumSize = pTrunc->sim()->hoz()->nL();
   }

   void AssocLegendreTransformOperatorTest::writeProjectors()
   {
      this->writeOperator("proj", &AssociatedLegendreOperator::proj, this->mProjExt);

      this->writeOperator("dProj", &AssociatedLegendreOperator::dProj, this->mProjExt);

      this->writeOperator("projS2Th", &AssociatedLegendreOperator::projS2Th, this->mProjExt);

      this->writeZOperator("projT2Th", &AssociatedLegendreOperator::projT2Th, this->mProjExt);

      this->writeZOperator("projS2Ph", &AssociatedLegendreOperator::projS2Ph, this->mProjExt);

      this->writeOperator("projT2Ph", &AssociatedLegendreOperator::projT2Ph, this->mProjExt);

      this->writeOperator("proj2GradTh", &AssociatedLegendreOperator::proj2GradTh, this->mProjExt);

      this->writeZOperator("proj2GradPh", &AssociatedLegendreOperator::proj2GradPh, this->mProjExt);
   }

   void AssocLegendreTransformOperatorTest::writeIntegrators()
   {
      this->writeOperator("intg", &AssociatedLegendreOperator::intg, this->mIntgExt);

      this->writeOperator("intgTh2S", &AssociatedLegendreOperator::intgTh2S, this->mIntgExt);

      this->writeZOperator("intgTh2T", &AssociatedLegendreOperator::intgTh2T, this->mIntgExt);

      this->writeZOperator("intgPh2S", &AssociatedLegendreOperator::intgPh2S, this->mIntgExt);

      this->writeOperator("intgPh2T", &AssociatedLegendreOperator::intgPh2T, this->mIntgExt);
   }
}
