/* \file POperator.hpp
 */

#ifndef POPERATOR_HPP
#define POPERATOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * @brief This class implements a polynomial operator expression of the type: P
    */
   class POperator
   {
      public:
         /// Typedef for the constant type
         typedef EPMFloat ConstantType;

         /// Typedef for full matrix type
         typedef Eigen::Matrix<ConstantType, Eigen::Dynamic, Eigen::Dynamic>  OutputMatrix;

         /**
          * @brief Constructor
          */
         POperator(const Matrix &poly);

         /**
          * @brief Destructor
          */
         virtual ~POperator();

         /**
          * @brief Get operator for the product (important if constant included)
          */
         const Matrix& productOp() const;

         /**
          * @brief Get operator
          */
         const Matrix& op() const;

         /**
          * @brief Multiplicative constant
          */
         EPMFloat c() const;
         
      protected:
         /**
          * @brief Storage reference base polynomial
          */
         const Matrix&  mrPoly;

      private:
   };

   inline const Matrix& POperator::productOp() const
   {
      return mrPoly;
   }

   inline const Matrix& POperator::op() const
   {
      return this->productOp();
   }

   inline EPMFloat POperator::c() const
   {
      return 1.0;
   }
}

#endif // POPERATOR_HPP
