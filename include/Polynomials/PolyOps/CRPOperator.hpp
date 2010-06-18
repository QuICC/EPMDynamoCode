/** \file CRPOperator.hpp
 */

#ifndef CRPOPERATOR_HPP
#define CRPOPERATOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * @brief This class implements a polynomial operator expression of the type: C P R
    *
    * \bug Needs cleaning and restructuring
    */
   template <typename TC> class CRPOperator
   {
      public:
         /// Typedef for the constant type
         typedef TC  ConstantType;

         /// Typedef for full matrix type
         typedef Eigen::Matrix<ConstantType, Eigen::Dynamic, Eigen::Dynamic>  OutputMatrix;

         /**
          * @brief Constructor
          */
         CRPOperator(const TC c, const Array &right, const Matrix &poly);

         /**
          * @brief Destructor
          */
         virtual ~CRPOperator();

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

         /**
          * @brief Storage reference right multipling array
          */
         const Array& mrR;

      private:
   };

   template <typename TC> CRPOperator<TC>::CRPOperator(const TC c, const Array &right, const Matrix &poly)
      : mC(c), mrPoly(poly), mrR(right)
   {
   }

   template <typename TC> CRPOperator<TC>::~CRPOperator()
   {
   }
   
   template <typename TC> inline typename CRPOperator<TC>::OutputMatrix CRPOperator<TC>::op() const
   {
      return mC*this->productOp();
   }
   
   template <typename TC> inline Matrix CRPOperator<TC>::productOp() const
   {
      return mrPoly*mrR.asDiagonal();
   }
   
   template <typename TC> inline TC CRPOperator<TC>::c() const
   {
      return mC;
   }

}

#endif // CRPOPERATOR_HPP
