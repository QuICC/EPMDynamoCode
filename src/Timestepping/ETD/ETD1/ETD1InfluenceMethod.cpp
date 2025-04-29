/** \file ETD1InfluenceMethod.cpp
 *  \brief Implementation of the ETD1 method (with influence matrix)
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Timestepping/ETD/ETD1/ETD1InfluenceKernel.hpp"

// Project includes
//

namespace EPMDynamo {

   ETD1InfluenceMethod::ETD1InfluenceMethod(EPMFloat a, EPMFloat b, const ETD1InfluenceMethod::BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0)
      : ETD1Method(a, b, basis, tsteps, pTrunc, hasL0), mInfluenceNBC(-2), mInfluence(pTrunc, basis, hasL0)
   {
      // Set the influence kernel
      SmartInfluenceKernel pIKernel(new ETDInfluenceKernel(1.0/this->mA, this->mETD1.pEtdF(1)));

      this->mpKernel = pIKernel;
   }

   void ETD1InfluenceMethod::addBC(SmartBC pBC)
   {
      // Propagate only the first boundary condition to standard implementation
      if(this->mInfluenceNBC == -2)
      {
         ETD1Method::addBC(pBC);
      }

      // Add boundary condition
      this->mInfluence.addBC(pBC);

      // Increment number of implemented influence matrix BCs
      ++this->mInfluenceNBC;
   }

   void ETD1InfluenceMethod::init()
   {
      // Initialise the timestep operators
      ETD1Method::init();

      // Initialise the influence matrix
      this->initInfluence();
   }

   void ETD1InfluenceMethod::updateTimeMatrices()
   {
      // Update the timestep operators
      ETD1Method::updateTimeMatrices();

      // Update the influence matrix solution
      this->updateInfluence();
   }

   void ETD1InfluenceMethod::doIteration(ETD1InfluenceMethod::ScalarType& rVar, ETD1InfluenceMethod::ScalarType& rNTerms)
   {
      // Solve influence matrix part
      this->mInfluence.solve(rNTerms);

      // Go on with normal timestep
      ETD1Method::doIteration(rVar, rNTerms);

      // Include kernel influence
      this->mInfluence.addKernel(rVar);
   }

   void ETD1InfluenceMethod::initInfluence()
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

   void ETD1InfluenceMethod::updateInfluence()
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
