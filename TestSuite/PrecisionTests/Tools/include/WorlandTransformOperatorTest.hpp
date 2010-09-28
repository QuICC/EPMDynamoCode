/** \file WorlandTransformOperatorTest.hpp
 *  \brief Worland specialisation of transform operators tests.
 */

#ifndef WORLANDTRANSFORMOPERATORTEST_HPP
#define WORLANDTRANSFORMOPERATORTEST_HPP

// System includes
//

// Project includes
//
#include "Polynomials/WorlandPolynomial.hpp"
#include "Polynomials/RadialBasis.hpp"

// TestSuite includes
//
#include "TransformOperatorTest.hpp"

namespace EPMDynamo {

typedef TorPolRadialOperator<WorlandPolynomial> WorlandOperator;

   /**
    * \brief Worland specialisation of transform operators tests.
    */
   class WorlandTransformOperatorTest: public TransformOperatorTest<RadialBasis<TorPolRadialOperator<WorlandPolynomial> > >
   {
      public:
         /**
         * @brief Constructs WorlandTransformOperatorTest.
         */
         WorlandTransformOperatorTest(SmartTruncation pTrunc);

         /**
         * @brief Simple empty destructor
         */
         virtual ~WorlandTransformOperatorTest() {};

         /**
          * @brief Write the projectors to file
          */
         void writeProjectors();
 
         /**
          * @brief Write the integrators to file
          */
         void writeIntegrators();        

         /**
          * @brief Get spectrum size depending on ordering index
          */
         int spectrumSize(const int index);

      protected:

      private:
   };

   inline int WorlandTransformOperatorTest::spectrumSize(const int index)
   {
      return this->mSpectrumSize;
   }

}

#endif // WORLANDTRANSFORMOPERATORTEST_HPP
