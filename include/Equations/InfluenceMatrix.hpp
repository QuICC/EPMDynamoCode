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
          * @brief Store kernel boundary values
          *
          * @param kernel Kernel influence to store
          * @param l Harmonic degree l
          */
         void storeKernelBC(const Array& kernel, const int l);

         /**
          * @brief Add kernel effect to given values
          *
          * @param rVar Solution to correct
          */
         void addKernel(ScalarType &rVar);

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

         /**
          * @brief Add boundary condition
          *
          * @param pBC Boundary condition
          */
         void addBC(SmartBC pBC);
         
      protected:

      private:
         /**
          * @brief Number of implemented BCs
          */
         int mOpBCs;

         /**
          * @brief Other boundary conditions
          */
         std::vector<SmartBC> mOtherBCs;

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
      : LaplacianBOperatorSet<TSimType, typename SimulationTraits<TSimType>::FactoredOpType>(basis, pTrunc), mOpBCs(-2)
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
      // Create the laplacian operators
      this->createOperators(-1.0);
   }

   template <typename TSimType> inline void InfluenceMatrix<TSimType>::addBC(SmartBC pBC)
   {
      // The first boundary conditions is implemented into the laplacian
      if(this->mOpBCs == -2)
      {
         LaplacianBOperatorSet<TSimType, typename SimulationTraits<TSimType>::FactoredOpType>::addBC(pBC);
      }
      // Other boundary conditions are stored in the influence matrix object
      else
      {
         this->mOtherBCs.push_back(pBC);
      }

      // increment number of implemented boundary conditions
      ++this->mOpBCs;
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

   template <typename TSimType> inline void InfluenceMatrix<TSimType>::storeKernelBC(const Array& kernel, const int l)
   {
      // Current implementation only works with a total 2 two BCs
      assert(this->mOtherBCs.size() == 1);

      // Compute the boundary value of the kernel
      EPMFloat bcVal = this->mOtherBCs.at(0)->getLHSBC(l).dot(kernel);

      // Protect against division by zero
      assert(bcVal != 0.0);

      // Store the rescaled kernel
      this->mSolutions.at(l) = kernel/bcVal;
   }

   template <typename TSimType> inline void InfluenceMatrix<TSimType>::addKernel(typename InfluenceMatrix<TSimType>::ScalarType &rVar)
   {
      // Get truncation information
      int nN = this->trunc()->sim()->rad()->nN();
      int nL = this->trunc()->local()->spec()->nL();
      const int l0 = rVar.minL();

      for(int l = l0; l < nL; ++l)
      {
         EPMComplex bcVal;

         for(int m =0; m <this->trunc()->local()->spec()->nM(l) ; ++m)
         {
            // Compute the boundary values of the timestep RHS (real and imaginary parts)
            bcVal.real() = this->mOtherBCs.at(0)->getLHSBC(l).dot(rVar.rLShell(l).col(m).real()); 
            bcVal.imag() = this->mOtherBCs.at(0)->getLHSBC(l).dot(rVar.rLShell(l).col(m).imag()); 

            // Loop over all radial coefficients
            for(int n=0; n < nN; ++n)
            {
               // Substract the kernel influence to give the right boundary condition
               rVar.rLShell(l).col(m)(n).real() -= bcVal.real()*this->mSolutions.at(l)(n);
               rVar.rLShell(l).col(m)(n).imag() -= bcVal.imag()*this->mSolutions.at(l)(n);
            }
         }
      }
   }

}

#endif // INFLUENCEMATRIX_HPP
