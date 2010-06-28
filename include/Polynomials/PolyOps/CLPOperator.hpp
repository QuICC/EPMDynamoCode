/** \file CLPOperator.hpp
 */

#ifndef CLPOPERATOR_HPP
#define CLPOPERATOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * @brief This class implements a polynomial operator expression of the type: C L P
    *
    * \epmBug Needs cleaning and restructuring
    */
   template <typename TC> class CLPOperator
   {
      public:
         /// Typedef for the constant type
         typedef TC  ConstantType;

         /// Typedef for full matrix type
         typedef Eigen::Matrix<ConstantType, Eigen::Dynamic, Eigen::Dynamic>  OutputMatrix;

         /**
          * @brief Constructor
          */
         CLPOperator(const TC c, const Array &left, const Matrix &poly);

         /**
          * @brief Destructor
          */
         virtual ~CLPOperator();

         /**
          * @brief Get operator for the product (important if constant included)
          */
         Matrix productOp() const;

         /**
          * @brief Get operator
          */
         OutputMatrix op() const;

         /**
          * @brief Multiplicative constant
          */
         TC c() const;
         
      protected:
         /**
          * @brief multiplicative constant
          */
         const TC mC;

         /**
          * @brief Storage reference base polynomial
          */
         const Matrix& mrPoly;

         /**
          * @brief Storage reference to the left multipling array
          */
         const Array& mrL;

      private:
   };

   template <typename TC> CLPOperator<TC>::CLPOperator(const TC c, const Array &left, const Matrix &poly)
      : mC(c), mrPoly(poly), mrL(left)
   {
   }

   template <typename TC> CLPOperator<TC>::~CLPOperator()
   {
   }
   
   template <typename TC> inline typename CLPOperator<TC>::OutputMatrix CLPOperator<TC>::op() const
   {
      return mC*this->productOp();
   }
   
   template <typename TC> inline Matrix CLPOperator<TC>::productOp() const
   {
      return mrL.asDiagonal()*mrPoly;
   }
   
   template <typename TC> inline TC CLPOperator<TC>::c() const
   {
      return mC;
   }

}

#endif // CLPOPERATOR_HPP
