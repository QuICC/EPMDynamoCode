/** \file CLRPOperator.hpp
 */

#ifndef CLRPOPERATOR_HPP
#define CLRPOPERATOR_HPP

// System includes
//

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * @brief This class implements a polynomial operator expression of the type: C L P R
    */
   template <typename TC> class CLRPOperator
   {
      public:
         /// Typedef for the constant type
         typedef TC  ConstantType;

         /// Typedef for full matrix type
         typedef Eigen::Matrix<ConstantType, Eigen::Dynamic, Eigen::Dynamic>  OutputMatrix;

         /**
          * @brief Constructor
          */
         CLRPOperator(const TC c, const Array &left, const Array &right, const Matrix &poly);

         /**
          * @brief Destructor
          */
         virtual ~CLRPOperator();

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
         const Array& mrL;

         /**
          * @brief Storage reference right multipling array
          */
         const Array& mrR;

      private:
   };

   template <typename TC> CLRPOperator<TC>::CLRPOperator(const TC c, const Array &left, const Array &right, const Matrix &poly)
      : mC(c), mrPoly(poly), mrL(left), mrR(right)
   {
   }

   template <typename TC> CLRPOperator<TC>::~CLRPOperator()
   {
   }
   
   template <typename TC> inline typename CLRPOperator<TC>::OutputMatrix CLRPOperator<TC>::op() const
   {
      return mC*this->productOp();
   }
   
   template <typename TC> inline Matrix CLRPOperator<TC>::productOp() const
   {
      return mrL.asDiagonal()*mrPoly*mrR.asDiagonal();
   }
   
   template <typename TC> inline TC CLRPOperator<TC>::c() const
   {
      return mC;
   }

}

#endif // CLRPOPERATOR_HPP
