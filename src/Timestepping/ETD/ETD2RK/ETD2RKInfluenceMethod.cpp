/** \file ETD2RKInfluenceMethod.cpp
 *  \brief Implementation of the ETD2RK method (with influence matrix)
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Timestepping/ETD/ETD2RK/ETD2RKInfluenceMethod.hpp"

// Project includes
//

namespace EPMDynamo {

   ETD2RKInfluenceMethod::ETD2RKInfluenceMethod(EPMFloat a, EPMFloat b, const ETD2RKInfluenceMethod::BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0)
      : ETD2RKMethod(a, b, basis, tsteps, pTrunc, hasL0), mInfluenceNBC(-2), mInfluence(pTrunc, basis, hasL0)
   {
      // Set the influence kernel
      SmartInfluenceKernel pIKernel(new ETDInfluenceKernel(1.0/this->mA, this->mETD2.pEtdF(1)));

      this->mpKernel = pIKernel;
   }

   void ETD2RKInfluenceMethod::addBC(SmartBC pBC)
   {
      // Propagate only the first boundary condition to standard implementation
      if(this->mInfluenceNBC == -2)
      {
         ETD2RKMethod::addBC(pBC);
      }

      // Add boundary condition
      this->mInfluence.addBC(pBC);

      // Increment number of implemented influence matrix BCs
      ++this->mInfluenceNBC;
   }

   void ETD2RKInfluenceMethod::init()
   {
      // Initialise the timestep operators
      ETD2RKMethod::init();

      // Initialise the influence matrix
      this->initInfluence();
   }

   void ETD2RKInfluenceMethod::updateTimeMatrices()
   {
      // Update the timestep operators
      ETD2RKMethod::updateTimeMatrices();

      // Update the influence matrix solution
      this->updateInfluence();
   }

   void ETD2RKInfluenceMethod::doIteration(ETD2RKInfluenceMethod::ScalarType& rVar, ETD2RKInfluenceMethod::ScalarType& rNTerms)
   {
      // Prepare values for iteration
      (*this->mCurrentIt)->prepare(rVar, rNTerms);

      // Solve influence matrix part
      this->mInfluence.solve(rNTerms);

      // Go on with normal timestep
      ETD2RKMethod::doIteration(rVar, rNTerms);

      // Include kernel influence
      this->mInfluence.addKernel(rVar);
   }

   void ETD2RKInfluenceMethod::initInfluence()
   {
      if(this->mInfluenceNBC  == 0)
      {
         // Initialise the operators
         this->mInfluence.initOperators();

         // Compute the operators
         this->mInfluence.computeOperators();
      } else
      {
         throw EPMException("ETD2RKInfluenceMethod::initInfluence", "Tried to initialise with wrong number of BCs");
      }
   }

   void ETD2RKInfluenceMethod::updateInfluence()
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
