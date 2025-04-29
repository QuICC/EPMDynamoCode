/** \file CPOperator.hpp
 */

#ifndef CPOPERATOR_HPP
#define CPOPERATOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * @brief This class implements a polynomial operator expression of the type: C P
    */
   template <typename TC> class CPOperator
   {
      public:
         /// Typedef for the constant type
         typedef TC  ConstantType;

         /// Typedef for full matrix type
         typedef Eigen::Matrix<ConstantType, Eigen::Dynamic, Eigen::Dynamic>  OutputMatrix;

         /**
          * @brief Constructor
          */
         CPOperator(const TC c, const Matrix &poly);

         /**
          * @brief Destructor
          */
         virtual ~CPOperator();

         /**
          * @brief Get operator
          */
         OutputMatrix op() const;

         /**
          * @brief Get operator for the product (important if constant included)
          */
         Matrix productOp() const;

         /**
          * @brief Multiplicative constant
          */
         TC c() const;

      protected:
         /**
          * @brief Multiplicative constant
          */
         const TC mC;

         /**
          * @brief Storage reference base polynomial
          */
         const Matrix& mrPoly;

      private:
   };

   template <typename TC> CPOperator<TC>::CPOperator(const TC c, const Matrix &poly)
      : mC(c), mrPoly(poly)
   {}

   template <typename TC> CPOperator<TC>::~CPOperator()
   {}

   template <typename TC> inline typename CPOperator<TC>::OutputMatrix CPOperator<TC>::op() const
   {
      return mC*this->productOp();
   }

   template <typename TC> inline Matrix CPOperator<TC>::productOp() const
   {
      return mrPoly;
   }

   template <typename TC> inline TC CPOperator<TC>::c() const
   {
      return mC;
   }

}

#endif // CPOPERATOR_HPP
