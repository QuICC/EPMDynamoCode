/** \file WorlandTransformOperatorTest.cpp
 *  \brief Definitions of the WorlandTransformOperatorTest methods.
 */

// System includes
//

// External includes
//

// Class include
//
#include "WorlandTransformOperatorTest.hpp"

// Project includes
//

namespace EPMDynamo {

   WorlandTransformOperatorTest::WorlandTransformOperatorTest(SmartTruncation pTrunc)
      : TransformOperatorTest<RadialBasis<TorPolRadialOperator<WorlandPolynomial> > >(pTrunc)
   {
      // Set prefix
      this->mPrefix = "wCheby_";

      // Set ordering string
      this->mOrderString = "_L";

      // Set max order index
      this->mOrderingN = pTrunc->sim()->hoz()->nL();

      // Set grid size
      this->mGridSize = pTrunc->sim()->rad()->nR();

      // Set spectrum size
      this->mSpectrumSize = pTrunc->sim()->rad()->nN();
   }

   void WorlandTransformOperatorTest::writeProjectors()
   {
      this->writeOperator("proj", &WorlandOperator::proj, this->mProjExt);

      this->writeOperator("dProj", &WorlandOperator::dProj, this->mProjExt);

      this->writeOperator("d2Proj", &WorlandOperator::d2Proj, this->mProjExt);

      this->writeOperator("proj2GradTP", &WorlandOperator::proj2GradTP, this->mProjExt);

      this->writeOperator("projPol2Q", &WorlandOperator::projPol2Q, this->mProjExt);

      this->writeOperator("projPol2S", &WorlandOperator::projPol2S, this->mProjExt);

      this->writeOperator("projTor2T", &WorlandOperator::projTor2T, this->mProjExt);

      this->writeOperator("projTor2CurlQ", &WorlandOperator::projTor2CurlQ, this->mProjExt);

      this->writeOperator("projTor2CurlS", &WorlandOperator::projTor2CurlS, this->mProjExt);

      this->writeOperator("projPol2CurlT", &WorlandOperator::projPol2CurlT, this->mProjExt);
   }

   void WorlandTransformOperatorTest::writeIntegrators()
   {
      this->writeOperator("intg", &WorlandOperator::intg, this->mIntgExt);

      this->writeOperator("intgQ2Pol", &WorlandOperator::intgQ2Pol, this->mIntgExt);

      this->writeOperator("intgT2Tor", &WorlandOperator::intgT2Tor, this->mIntgExt);

      this->writeOperator("intgQ2CurlCurlProj", &WorlandOperator::intgQ2CurlCurlProj, this->mIntgExt);

      this->writeOperator("intgS2CurlCurlProj", &WorlandOperator::intgS2CurlCurlProj, this->mIntgExt);

      this->writeOperator("intgT2CurlProj", &WorlandOperator::intgT2CurlProj, this->mIntgExt);
   }
}
