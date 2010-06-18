/** \file QuadratureMultiplier.hpp
 *  \brief Implementation of a general quadrature transform
 */

#ifndef QUADRATUREMULTIPLIER_HPP
#define QUADRATUREMULTIPLIER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "GeneralScalars/FDSHShell.hpp"
#include "Transforms/Polynomial/QuadratureBasis.hpp"
#include "Polynomials/PolynomialOperator.hpp"

namespace EPMDynamo {
   
   /**
    * @brief This class contains the methods that are required for quadrature transforms
    *
    * The class provides the methods required for the quadrature integration and
    * quadrature summation which are required for the transforms.
    *
    * \tparam TBasis Basis of the quadrature
    *
    * \bug Review this and simplify if possible. Changes will depend of new implementation of Polynomial Operators.
    */
   template<typename TBasis> class QuadratureMultiplier: public QuadratureBasis<TBasis>
   {
      public:
         /// Typedef to get access to the type of the underlying Basis
         typedef TBasis    PolynomialBasis;

         /**
         * @brief Constructs the inherited basis
         *
         * This class only contains methods, the basis object gets inherited
         *
         * @param params   Array of the fixed parameters
         * @param nX       Size of the grid
         * @param nN       Size of the polynomial
         */
         QuadratureMultiplier(const ArrayI& params, const int nX, const int nN);

         /**
         * @brief Simple empty destructor
         */
         virtual ~QuadratureMultiplier() {};

         /**
          * @brief Fully templated multiplication by polynomial
          */
         template<class TOp, class T, typename T1, typename T2, typename TExp> void multL(T1 &rOutVal,const T2 &inVal, const PolynomialOperator<TExp>& (T::*Tpr)() const );
         /**
          * @brief Overloaded multiplication by polynomial for FDSH transforms (forward)
          */
         template<class TOp, class T, typename TExp> void multL(FDSHShell &rOutVal,const MatrixZ &inVal, const PolynomialOperator<TExp>& (T::*Tpr)() const );
         /**
          * @brief Overloaded multiplication by polynomial for FDSH transforms (backward)
          */
         template<class TOp, class T, typename TExp> void multL(MatrixZ &rOutVal,const FDSHShell &inVal, const PolynomialOperator<TExp>& (T::*Tpr)() const );

         /**
          * @brief Fully templated multiplication by polynomial
          */
         template<class TOp, class T, typename T1, typename T2, typename TExp> void multM(T1 &rOutVal,const T2 &inVal, const PolynomialOperator<TExp>& (T::*Tpr)() const );

      protected:

      private:
   };

   template<typename TBasis>  QuadratureMultiplier<TBasis>::QuadratureMultiplier(const ArrayI& params, const int nX, const int nN)
      : QuadratureBasis<TBasis>(params, nX, nN)
   {
   }

   template<typename TBasis> template<class TOp, class T, typename T1, typename T2, typename TExp> inline void  QuadratureMultiplier<TBasis>::multL(T1 &rOutVal, const T2 &inVal, const PolynomialOperator<TExp>& (T::*Tptr)() const)
   {
      // Perform integration over polynomials
      for(int k=0; k < this->mN; ++k)
      {
         TOp::eval(rOutVal.rLShell(k) , (this->mBasis.at(k).*Tptr)().productOp().transpose() * inVal.lshell(k)*(this->mBasis.at(k).*Tptr)().c());
      }
   }

   template<typename TBasis> template<class TOp, class T, typename TExp> inline void  QuadratureMultiplier<TBasis>::multL(FDSHShell &rOutVal, const MatrixZ &inVal, const PolynomialOperator<TExp>& (T::*Tptr)() const)
   {
      // Perform integration over polynomials
      for(int k=0; k < this->mN; ++k)
      {
         TOp::eval(rOutVal.rOrder(k), (this->mBasis.at(k).*Tptr)().productOp().transpose() * inVal.col(k)*(this->mBasis.at(k).*Tptr)().c());
      }
   }

   template<typename TBasis> template<class TOp, class T, typename TExp> inline void  QuadratureMultiplier<TBasis>::multL(MatrixZ &rOutVal, const FDSHShell &inVal, const PolynomialOperator<TExp>& (T::*Tptr)() const)
   {
      // Perform integration over polynomials
      for(int k=0; k < this->mN; ++k)
      {
         TOp::eval(rOutVal, k, (this->mBasis.at(k).*Tptr)().productOp().transpose() * inVal.order(k)*(this->mBasis.at(k).*Tptr)().c());
      }
   }

   template<typename TBasis> template<class TOp, class T, typename T1, typename T2, typename TExp> inline void  QuadratureMultiplier<TBasis>::multM(T1 &rOutVal, const T2 &inVal, const PolynomialOperator<TExp>& (T::*Tptr)() const)
   {
      // Perform integration over polynomials
      for(int k=0; k < this->mN; ++k)
      {
         TOp::eval(rOutVal.rMShell(k) , (this->mBasis.at(k).*Tptr)().productOp().transpose() * inVal.mshell(k)*(this->mBasis.at(k).*Tptr)().c());
      }
   }
}

#endif // QUADRATUREMULTIPLIER_HPP
