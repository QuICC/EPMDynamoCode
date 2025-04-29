/** \file AssocLegendreTransformOperatorTest.hpp
 *  \brief Associated Legendre specialisation of transform operators tests.
 */

#ifndef ASSOCLEGENDRETRANSFORMOPERATORTEST_HPP
#define ASSOCLEGENDRETRANSFORMOPERATORTEST_HPP

// System includes
//

// Project includes
//
#include "Polynomials/AssocLegendreBasis.hpp"

// TestSuite includes
//
#include "TransformOperatorTest.hpp"

namespace EPMDynamo {

   /**
    * \brief Worland specialisation of transform operators tests.
    */
   class AssocLegendreTransformOperatorTest: public TransformOperatorTest<AssocLegendreBasis>
   {
      public:
         /**
         * @brief Constructs AssocLegendreTransformOperatorTest.
         */
         AssocLegendreTransformOperatorTest(SmartTruncation pTrunc);

         /**
         * @brief Simple empty destructor
         */
         virtual ~AssocLegendreTransformOperatorTest() {};

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

   inline int AssocLegendreTransformOperatorTest::spectrumSize(const int index)
   {
      return this->mSpectrumSize - index;
   }

}

#endif // ASSOCLEGENDRETRANSFORMOPERATORTEST_HPP
