/** \file RPOperator.hpp
 */

#ifndef RPOPERATOR_HPP
#define RPOPERATOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * @brief This class implements a polynomial operator expression of the type: P R
    */
   class RPOperator
   {
      public:
         /// Typedef for the constant type
         typedef EPMFloat ConstantType;

         /// Typedef for full matrix type
         typedef Eigen::Matrix<ConstantType, Eigen::Dynamic, Eigen::Dynamic>  OutputMatrix;

         /**
          * @brief Constructor
          */
         RPOperator(const Array &right, const Matrix &poly);

         /**
          * @brief Destructor
          */
         virtual ~RPOperator();

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
          * @brief Storage reference right multipling array
          */
         const Array& mrR;

      private:
   };

   inline Matrix RPOperator::productOp() const
   {
      return (mrPoly*mrR.asDiagonal());
   }

   inline Matrix RPOperator::op() const
   {
      return this->productOp();
   }

   inline EPMFloat RPOperator::c() const
   {
      return 1.0;
   }

}

#endif // RPOPERATOR_HPP
