/** \file LPOperator.hpp
 */

#ifndef LPOPERATOR_HPP
#define LPOPERATOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * @brief This class implements a polynomial operator expression of the type: L P
    */
   class LPOperator
   {
      public:
         /// Typedef for the constant type
         typedef EPMFloat ConstantType;

         /// Typedef for full matrix type
         typedef Eigen::Matrix<ConstantType, Eigen::Dynamic, Eigen::Dynamic>  OutputMatrix;

         /**
          * @brief Constructor
          */
         LPOperator(const Array &left, const Matrix &poly);

         /**
          * @brief Destructor
          */
         virtual ~LPOperator();

         /**
          * @brief Get operator for the product (important if constant included)
          */
         Matrix productOp() const;

         /**
          * @brief Get operator
          */
         Matrix op() const;

         /**
          * @brief Multiplicative constant
          */
         EPMFloat c() const;
         
      protected:
         /**
          * @brief Storage reference base polynomial
          */
         const Matrix& mrPoly;

         /**
          * @brief Storage reference left multipling array
          */
         const Array& mrL;

      private:
   };

   inline Matrix LPOperator::productOp() const
   {
      return mrL.asDiagonal()*mrPoly;
   }

   inline Matrix LPOperator::op() const
   {
      return this->productOp();
   }

   inline EPMFloat LPOperator::c() const
   {
      return 1.0;
   }

}

#endif // LPOPERATOR_HPP
