/** \file AssociatedLegendreOperator.hpp
 *  \brief Implementation of a set of operators on the associated Legendre polynomial
 */

#ifndef ASSOCIATEDLEGENDREOPERATOR_HPP
#define ASSOCIATEDLEGENDREOPERATOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Polynomials/AssociatedLegendreOperatorBase.hpp"
#include "Polynomials/PolynomialOperator.hpp"
#include "Polynomials/PolyOps/PolyOpsTypes.hpp"
#include "General/MathConstants.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a set of operators on the associated Legendre polynomial
    *
    * \bug Review this to remove the need of the complicated copy/assignement operators
    * \bug Rewrite the polynomial operators in a cleaner and simpler way
    */
   class AssociatedLegendreOperator: public AssociatedLegendreOperatorBase
   {
      public:

         /**
          * @brief Constructs the polynomial and the operators from scratch.
          *
          * This constructors computes the values of the grid and the weights and then
          * computes the polynomial values based on these
          *
          * @param nX      Size of the grid
          * @param nL      Number of l modes
          * @param pM      Polynomial m degree
          */
         AssociatedLegendreOperator(const int pM, const int nX, const int nL);

         /**
          * @brief Constructs the polynomial and the operators based on the given grid points and weights
          *
          * The values of the polynomial are computed on the given grid and weighted with
          * the given weights.
          *
          * @param grid    Provided grid array
          * @param nL      Number of l modes
          * @param weights Provided weights array
          * @param pM      Polynomial m degree
          */
         AssociatedLegendreOperator(const int pM, const SmartArray grid, const int nL, const SmartArray weights);

         /**
          * @brief Need copy constructor to avoid shallow copy!
          *
          * @param poly Polynomial to copy
          */
         AssociatedLegendreOperator(const AssociatedLegendreOperator &poly);

         /**
          * @brief Need copy assignement to avoid shallow copy!
          *
          * @param poly Polynomial to assign
          */
         AssociatedLegendreOperator&  operator=(const AssociatedLegendreOperator &poly);

         /**
          * @brief Simple empty destructor
          */
         virtual ~AssociatedLegendreOperator() {};

         //
         // Integration operators
         //

         /**
          * @brief Get the integration operator
          */
         const PolynomialOperator<POperator>& intg() const;

         /**
          * @brief Get the integration operator from theta component to the QST S 
          *          component
          */
         const PolynomialOperator<RPOperator>& intgTh2S() const;

         /**
          * @brief Get the integration operator from phi component to the QST S 
          *          component
          */
         const PolynomialOperator<CLRPOperator<EPMComplex> >& intgPh2S() const;

         /**
          * @brief Get the integration operator from theta component to the QST T 
          *          component
          */
         const PolynomialOperator<CLRPOperator<EPMComplex> >& intgTh2T() const;

         /**
          * @brief Get the integration operator from phi component to the QST T 
          *          component
          */
         const PolynomialOperator<RPOperator>& intgPh2T() const;

         //
         // Projection operators
         //

         /**
          * @brief Get the projection operator expansion
          */
         const PolynomialOperator<POperator>& proj() const;

         /**
          * @brief Get the projection operator of the first derivative 
          */
         const PolynomialOperator<POperator>& dProj() const;

         /**
          * @brief Get the projection operator from the QST S component to the theta
          *          component 
          */
         const PolynomialOperator<LPOperator>& projS2Th() const;

         /**
          * @brief Get the projection operator from the QST S component to the phi
          *          component 
          */
         const PolynomialOperator<CLRPOperator<EPMComplex> >& projS2Ph() const;

         /**
          * @brief Get the projection operator from the QST T component to the theta
          *          component 
          */
         const PolynomialOperator<CLRPOperator<EPMComplex> >& projT2Th() const;

         /**
          * @brief Get the projection operator from the QST T component to the phi
          *          component 
          */
         const PolynomialOperator<LPOperator>& projT2Ph() const;

         /**
          * @brief Get the projection operator to the theta component  of the gradient
          */
         const PolynomialOperator<POperator>& proj2GradTh() const;

         /**
          * @brief Get the projection operator to the phi component  of the gradient
          */
         const PolynomialOperator<CRPOperator<EPMComplex> >& proj2GradPh() const;
         
      protected:

      private:
         /**
          * @brief Stored smart ZRP Operator
          */
         SmartZRPOp  mpProj2GradPh;

         /**
          * @brief Stored smart P Operator
          */
         SmartPOp mpProj2GradTh;

         /**
          * @brief Stored smart LP Operator
          */
         SmartLPOp   mpProjT2Ph;

         /**
          * @brief Stored smart ZLRP Operator
          */
         SmartZLRPOp   mpProjT2Th;

         /**
          * @brief Stored smart ZLRP Operator
          */
         SmartZLRPOp   mpProjS2Ph;

         /**
          * @brief Stored smart LP Operator
          */
         SmartLPOp   mpProjS2Th;

         /**
          * @brief Stored smart P Operator
          */
         SmartPOp    mpDProj;

         /**
          * @brief Stored smart P Operator
          */
         SmartPOp    mpProj;

         /**
          * @brief Stored smart RP Operator
          */
         SmartRPOp      mpIntgPh2T;

         /**
          * @brief Stored smart ZLRP Operator
          */
         SmartZLRPOp  mpIntgTh2T;

         /**
          * @brief Stored smart ZLRP Operator
          */
         SmartZLRPOp   mpIntgPh2S;

         /**
          * @brief Stored smart RP Operator
          */
         SmartRPOp  mpIntgTh2S;

         /**
          * @brief Stored smart ZP Operator
          */
         SmartPOp  mpIntg;
         
         /**
          * @brief Init all the PolynomialOperator<*> objects
          */
         void initPOperators();

   };

   //
   // Inline definitions of the integrators
   //

   inline const PolynomialOperator<POperator>& AssociatedLegendreOperator::intg() const
   {
      return *this->mpIntg;
   }

   inline const PolynomialOperator<RPOperator>& AssociatedLegendreOperator::intgTh2S() const
   {
      return *this->mpIntgTh2S;
   }

   inline const PolynomialOperator<CLRPOperator<EPMComplex> >& AssociatedLegendreOperator::intgPh2S() const
   {
      return *this->mpIntgPh2S;
   }

   inline const PolynomialOperator<CLRPOperator<EPMComplex> >& AssociatedLegendreOperator::intgTh2T() const
   {
      return *this->mpIntgTh2T;
   }

   inline const PolynomialOperator<RPOperator>& AssociatedLegendreOperator::intgPh2T() const
   {
      return *this->mpIntgPh2T;
   }

   //
   // Inline definitions of the projectors
   //

   inline const PolynomialOperator<POperator>& AssociatedLegendreOperator::proj() const
   {
      return *this->mpProj;
   }

   inline const PolynomialOperator<POperator>& AssociatedLegendreOperator::dProj() const
   {
      return *this->mpDProj;
   }

   inline const PolynomialOperator<LPOperator>& AssociatedLegendreOperator::projS2Th() const
   {
      return *this->mpProjS2Th;
   }

   inline const PolynomialOperator<CLRPOperator<EPMComplex> >& AssociatedLegendreOperator::projS2Ph() const
   {
      return *this->mpProjS2Ph;
   }

   inline const PolynomialOperator<CLRPOperator<EPMComplex> >& AssociatedLegendreOperator::projT2Th() const
   {
      return *this->mpProjT2Th;
   }

   inline const PolynomialOperator<LPOperator>& AssociatedLegendreOperator::projT2Ph() const
   {
      return *this->mpProjT2Ph;
   }

   inline const PolynomialOperator<POperator>& AssociatedLegendreOperator::proj2GradTh() const
   {
      return *this->mpProj2GradTh;
   }

   inline const PolynomialOperator<CRPOperator<EPMComplex> >& AssociatedLegendreOperator::proj2GradPh() const
   {
      return *this->mpProj2GradPh;
   }

}

#endif // ASSOCIATEDLEGENDREOPERATOR_HPP
