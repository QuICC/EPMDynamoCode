/** \file LRPOperator.hpp
 */

#ifndef LRPOPERATOR_HPP
#define LRPOPERATOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    *
    * \epmBug Needs cleaning and restructuring
    */
   class LRPOperator
   {
      public:
         /// Typedef for the constant type
         typedef EPMFloat ConstantType;

         /// Typedef for full matrix type
         typedef Eigen::Matrix<ConstantType, Eigen::Dynamic, Eigen::Dynamic>  OutputMatrix;

         /**
          * @brief Constructor
          */
         LRPOperator(const Array &left, const Array &right, const Matrix &poly);

         /**
          * @brief Destructor
          */
         virtual ~LRPOperator();

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
          * @brief Storage reference left multipling array
          */
         const Array& mrL;

         /**
          * @brief Storage reference right multipling array
          */
         const Array& mrR;

         /**
          * @brief Storage reference base polynomial
          */
         const Matrix& mrPoly;

      private:
   };

   inline Matrix LRPOperator::productOp() const
   {
      return mrL.asDiagonal()*mrPoly*mrR.asDiagonal();
   }

   inline Matrix LRPOperator::op() const
   {
      return this->productOp();
   }

   inline EPMFloat LRPOperator::c() const
   {
      return 1.0;
   }

}

#endif // LRPOPERATOR_HPP
