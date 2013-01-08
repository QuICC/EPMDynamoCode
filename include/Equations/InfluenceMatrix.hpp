/** \file InfluenceMatrix.hpp
 *  \brief Implementation of the influence matrix approach
 */

#ifndef INFLUENCEMATRIX_HPP
#define INFLUENCEMATRIX_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//
#include <vector>

// External includes
//
#include "LAPACK_Iface.hpp"

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"
#include "Operators/LaplacianBOperatorSet.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the influence matrix approach for the 4th order Navier-Stokes equation
    */
   class InfluenceMatrix: public LaplacianBOperatorSet<SimulationConfig::FactoredOpType>
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef SimulationConfig::NumericalScheme::RadialBasisType    BasisType;

         /// Typedef from Simulation trait to local truncation type
         typedef SimulationConfig::NumericalScheme::ScalarType    ScalarType;

         /**
          * @brief Constructor
          *
          * @param pTrunc Truncation information
          * @param basis Reference to the basis used for the operators
          * @param hasL0 Is the l=0 mode required?
          */
         InfluenceMatrix(SmartTruncation pTrunc, const BasisType &basis, bool hasL0);

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

   inline void InfluenceMatrix::addBC(SmartBC pBC)
   {
      // The first boundary conditions is implemented into the laplacian
      if(this->mOpBCs == -2)
      {
         LaplacianBOperatorSet<SimulationConfig::FactoredOpType>::addBC(pBC);
      }
      // Other boundary conditions are stored in the influence matrix object
      else
      {
         this->mOtherBCs.push_back(pBC);
      }

      // increment number of implemented boundary conditions
      ++this->mOpBCs;
   }

   inline void InfluenceMatrix::solve(InfluenceMatrix::ScalarType &rVar)
   {
      int degrees = this->trunc()->local()->spec()->nL();
      const int l0 = rVar.minL();

      for(int l = l0; l < degrees; ++l)
      {
         this->solveOrders(rVar.rLShell(l), l);
      }
   }

   inline void InfluenceMatrix::storeKernelBC(const Array& kernel, const int l)
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

   inline void InfluenceMatrix::addKernel(InfluenceMatrix::ScalarType &rVar)
   {
      // Get truncation information
      int nN = this->trunc()->sim()->rad()->nN();
      int nL = this->trunc()->local()->spec()->nL();
      const int l0 = rVar.minL();

      int l_;
      int m_;

      for(int l = l0; l < nL; ++l)
      {
         EPMComplex bcVal;

         l_ = this->trunc()->local()->spec()->lArray()(l);

         for(int m =0; m <this->trunc()->local()->spec()->nM(l) ; ++m)
         {
            m_ = this->trunc()->local()->spec()->mArray(l)(m);

            // Compute the boundary values of the timestep RHS (real and imaginary parts)
            bcVal.real() = this->mOtherBCs.at(0)->getLHSBC(l).dot(rVar.rLShell(l).col(m).real()); 
            bcVal.imag() = this->mOtherBCs.at(0)->getLHSBC(l).dot(rVar.rLShell(l).col(m).imag());

            // Loop over all radial coefficients
            for(int n=0; n < nN; ++n)
            {
               // Substract the kernel influence to give the right boundary condition
               rVar.rLShell(l).col(m)(n).real() += (this->mOtherBCs.at(0)->getRealRHSBC(l_,m_) - bcVal.real())*this->mSolutions.at(l)(n);
               rVar.rLShell(l).col(m)(n).imag() += (this->mOtherBCs.at(0)->getImagRHSBC(l_,m_) - bcVal.imag())*this->mSolutions.at(l)(n);
            }
         }
      }
   }

}

#endif // INFLUENCEMATRIX_HPP
