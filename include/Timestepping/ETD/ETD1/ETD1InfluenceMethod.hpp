/** \file ETD1InfluenceMethod.hpp
 *  \brief Implementation of the ETD1 method (with influence matrix)
 */

#ifndef ETD1INFLUENCEMETHOD_HPP
#define ETD1INFLUENCEMETHOD_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "Timestepping/TimestepParameters.hpp"
#include "Timestepping/ETD/ETD1/ETD1Method.hpp"
#include "Timestepping/ETD/ETDInfluenceKernel.hpp"
#include "Config/SimulationConfig.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the theta method (without influence matrix)
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ETD1InfluenceMethod: public ETD1Method<TSimType>
   {
      public:
         /// Typedef from Simulation trait to local radial basis type
         typedef typename TSimType::RadialBasisType    BasisType;

         /// Typedef from Simulation trait to local scalar type
         typedef typename TSimType::ScalarType    ScalarType;

         /// Typedef for a smart ETD1 influence kernel
         typedef EPMSHARED_PTR<ETDInfluenceKernel<TSimType> > SmartInfluenceKernel;

         /**
          * @brief Constructor
          *
          * @param a Coefficient \f$a\f$ of timestep scheme
          * @param b Coefficient \f$b\f$ of timestep scheme
          * @param basis Reference to the basis used for the operators
          * @param tsteps Timestep parameters
          * @param pTrunc Truncation information
          * @param hasL0 Is the l=0 mode required?
          */
         ETD1InfluenceMethod(EPMFloat a, EPMFloat b, const BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0);

         /**
          * @brief Destructor
          */
         virtual ~ETD1InfluenceMethod() {};

         /**
          * @brief Add boundary condition
          *
          * Only the first boundary condition will be propagated to the standard timestep structure. But all boundray conditions
          * are inserted into the influence matrix.
          *
          * @param pBC Boundary condition
          */
         void addBC(SmartBC pBC);

         /**
          * @brief Initialise the ETD1 influence matrix method
          */
         void init();
         
      protected:

         /**
          * @brief Update the timestep matrices after a timestep change
          */
         void updateTimeMatrices();

         /**
          * @brief Compute next scheme iteration
          *
          * @param rVar Input/Output variable
          * @param rNTerms New non linear terms
          */
         void doIteration(ScalarType& rVar, ScalarType& rNTerms);

      private:
         /**
          * @brief Number of boundary conditions required for influence matrix
          */
         int mInfluenceNBC;

         /**
          * @brief Influence matrix
          */
         InfluenceMatrix<TSimType> mInfluence;

         /**
          * @brief Computation of the influence kernel
          */
         SmartInfluenceKernel  mpKernel;

         /**
          * @brief Initialise the influenc matrix
          */
         void initInfluence();

         /**
          * @brief update the influence matrix solution
          */
         void updateInfluence();
   };

   template <typename TSimType> ETD1InfluenceMethod<TSimType>::ETD1InfluenceMethod(EPMFloat a, EPMFloat b, const typename ETD1InfluenceMethod<TSimType>::BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0)
      : ETD1Method<TSimType>(a, b, basis, tsteps, pTrunc, hasL0), mInfluenceNBC(-2), mInfluence(pTrunc, basis, hasL0)
   {
      // Set the influence kernel
      SmartInfluenceKernel pIKernel(new ETDInfluenceKernel<TSimType>(1.0/this->mA, this->mETD1.pEtdF(1)));

      this->mpKernel = pIKernel;
   }

   template <typename TSimType> void ETD1InfluenceMethod<TSimType>::addBC(SmartBC pBC)
   {
      // Propagate only the first boundary condition to standard implementation
      if(this->mInfluenceNBC == -2)
      {
         ETD1Method<TSimType>::addBC(pBC);
      }

      // Add boundary condition
      this->mInfluence.addBC(pBC);

      // Increment number of implemented influence matrix BCs
      ++this->mInfluenceNBC;
   }

   template <typename TSimType> void ETD1InfluenceMethod<TSimType>::init()
   {
      // Initialise the timestep operators
      ETD1Method<TSimType>::init();

      // Initialise the influence matrix
      this->initInfluence();
   }

   template <typename TSimType> void ETD1InfluenceMethod<TSimType>::updateTimeMatrices()
   {
      // Update the timestep operators
      ETD1Method<TSimType>::updateTimeMatrices();

      // Update the influence matrix solution
      this->updateInfluence();
   }

   template <typename TSimType> void ETD1InfluenceMethod<TSimType>::doIteration(typename ETD1InfluenceMethod<TSimType>::ScalarType& rVar, typename ETD1InfluenceMethod<TSimType>::ScalarType& rNTerms)
   {
      // Solve influence matrix part
      this->mInfluence.solve(rNTerms);

      // Go on with normal timestep
      ETD1Method<TSimType>::doIteration(rVar, rNTerms);

      // Include kernel influence
      this->mInfluence.addKernel(rVar);
   }

   template <typename TSimType> void ETD1InfluenceMethod<TSimType>::initInfluence()
   {
      if(this->mInfluenceNBC  == 0)
      {
         // Initialise the operators
         this->mInfluence.initOperators();

         // Compute the operators
         this->mInfluence.computeOperators();
      } else
      {
         throw EPMException("ETD1InfluenceMethod::initInfluence", "Tried to initialise with wrong number of BCs");
      }
   }

   template <typename TSimType> void ETD1InfluenceMethod<TSimType>::updateInfluence()
   {
      // Get size of radial truncation
      int nN = this->mInfluence.trunc()->sim()->rad()->nN();
      // Get number of harmonic degrees
      int nL = this->mInfluence.trunc()->local()->spec()->nL();

      // Get minimal degree index (not l=0)
      int l0 = this->mInfluence.minL();

      // Create temporary storage
      Array tmp(nN);

      // loop over degrees 
      for(int l = l0; l < nL; ++l)
      {
         // Initialise influence matrix solution to kernel r^l
         tmp.setConstant(0.0);
         tmp(0) = 1.0;

         // Compute the kernel influence
         this->mpKernel->computeInfluence(tmp, l);

         // Store solution from influence matrix
         this->mInfluence.storeKernelBC(tmp, l);
      }
   }

}

#endif // ETD1INFLUENCEMETHOD_HPP
