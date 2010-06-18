/** \file InfluenceMatrix.hpp
 *  \brief Implementation of the influence matrix approach
 */

#ifndef INFLUENCEMATRIX_HPP
#define INFLUENCEMATRIX_HPP

// System includes
//
#include <vector>

// External includes
//
#include "LAPACK_Iface.hpp"

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Simulations/Traits/SimulationTraits.hpp"
#include "Domain/Truncation.hpp"
#include "Operators/LaplacianBOperatorSet.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the influence matrix approach for the 4th order Navier-Stokes equation
    *
    * \tparam TSimType Type of the simulation
    *
    * \bug Review computations
    */
   template <typename TSimType> class InfluenceMatrix: public LaplacianBOperatorSet<TSimType, typename SimulationTraits<TSimType>::FactoredOpType>
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::RadialBasisType    BasisType;

         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::ScalarType    ScalarType;

         /**
          * @brief Constructor
          *
          * @param pTrunc Truncation information
          * @param basis Reference to the basis used for the operators
          */
         InfluenceMatrix(SmartTruncation pTrunc, const BasisType &basis);

         /**
          * @brief Simple empty destructor
          */
         virtual ~InfluenceMatrix() {};

         /**
          * @brief Store solution
          *
          * @param sol Solution to store
          * @param l Harmonic degree l
          *
          * \bug IS NOT BOUNDARY CONDITION AWARE FOR THE MOMENT
          */
         void storeSolution(const Array& sol, const int l);

         /**
          * @brief Correct timestep solution
          *
          * @param rVar Intermediate solution to correcte
          *
          * \bug IS NOT BOUNDARY CONDITION AWARE FOR THE MOMENT
          */
         void correctSolution(ScalarType &rVar);

         /**
          * @brief Solve influence matrix equation for unkown
          *
          * @param rVar Storage for solution of linear equation
          */
         void solve(ScalarType &rVar);

         /**
          * @brief Compute the operators
          *
          * Computation is done in three steps.
          *    1) Fill in laplacian
          *    2) Compute LU factorisation
          */
         void computeOperators();
         
      protected:

      private:
         /**
          * @brief Vector of solution arrays
          */
         std::vector<Array>   mSolutions;

         /**
          * @brief Initialise the solution vector
          */
         void initSolutions();
   };

   template <typename TSimType> InfluenceMatrix<TSimType>::InfluenceMatrix(SmartTruncation pTrunc, const typename InfluenceMatrix<TSimType>::BasisType &basis)
      : LaplacianBOperatorSet<TSimType, typename SimulationTraits<TSimType>::FactoredOpType>(basis, pTrunc)
   {
      this->initSolutions();
   }

   template <typename TSimType> void InfluenceMatrix<TSimType>::initSolutions()
   {
      int nN = this->trunc()->sim()->rad()->nN();
      int nL = this->trunc()->local()->spec()->nL();

      for(int l = 0; l < nL; ++l)
      {
         this->mSolutions.push_back(Array(nN));
         this->mSolutions.at(l).setConstant(0.0);
      }
   }

   template <typename TSimType> void InfluenceMatrix<TSimType>::computeOperators()
   {
      this->createOperators(-1.0);
   }

   template <typename TSimType> inline void InfluenceMatrix<TSimType>::solve(typename InfluenceMatrix<TSimType>::ScalarType &rVar)
   {
      int degrees = this->trunc()->local()->spec()->nL();
      const int l0 = rVar.minL();

      for(int l = l0; l < degrees; ++l)
      {
         this->solveOrders(rVar.rLShell(l), l);
      }
   }

   template <typename TSimType> inline void InfluenceMatrix<TSimType>::storeSolution(const Array& sol, const int l)
   {
      EPMFloat bcVal = this->mrBasis.at(l).bpoly().dot(sol);

      this->mSolutions.at(l) = sol/bcVal;
   }

   template <typename TSimType> inline void InfluenceMatrix<TSimType>::correctSolution(typename InfluenceMatrix<TSimType>::ScalarType &rVar)
   {
      int nN = this->trunc()->sim()->rad()->nN();
      int nL = this->trunc()->local()->spec()->nL();
      const int l0 = rVar.minL();

      for(int l = l0; l < nL; ++l)
      {
         EPMComplex bcVal;

         for(int m =0; m <this->trunc()->local()->spec()->nM(l) ; ++m)
         {
            bcVal.real() = this->mrBasis.at(l).bpoly().dot(rVar.rLShell(l).col(m).real()); 
            bcVal.imag() = this->mrBasis.at(l).bpoly().dot(rVar.rLShell(l).col(m).imag()); 
            for(int n=0; n < nN; ++n)
            {
               rVar.rLShell(l).col(m)(n).real() -= bcVal.real()*mSolutions.at(l)(n);
               rVar.rLShell(l).col(m)(n).imag() -= bcVal.imag()*mSolutions.at(l)(n);
            }
         }
      }
   }

}

#endif // INFLUENCEMATRIX_HPP
