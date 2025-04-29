/** \file IterativeSchemeBase.hpp
 *  \brief Base of the implementation of an iterative timestep scheme
 */

#ifndef ITERATIVESCHEMEBASE_HPP
#define ITERATIVESCHEMEBASE_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"
#include "Config/NumericalSchemeInc.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "BoundaryConditions/BoundaryCondition.hpp"
#include "Timestepping/TimestepParameters.hpp"
#include "Timestepping/TimestepSchemeBase.hpp"
#include "Timestepping/SchemeIteration.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the implementation of an iterative timestep scheme
    */
   class IterativeSchemeBase: public TimestepSchemeBase
   {
      public:
         /// Typedef from Simulation trait to local radial basis type
         typedef SimulationConfig::NumericalScheme::RadialBasisType    BasisType;

         /// Typedef from Simulation trait to local scalar type
         typedef SimulationConfig::NumericalScheme::ScalarType    ScalarType;

         /// Typedef for a smart pointer to a scalar type
         typedef EPMSHARED_PTR<ScalarType>   SmartScalarType;

         /// Typedef for the smart pointer to an general iteration
         typedef EPMSHARED_PTR<SchemeIteration >   SmartIteration;

         /**
          * @brief Constructor
          *
          * @param a Coefficient \f$a\f$ of timestep scheme
          * @param b Coefficient \f$b\f$ of timestep scheme
          * @param basis Reference to the basis used for the operators
          * @param tsteps Timestep parameters
          * @param pTrunc Truncation information
          * @param hasL0 Is l=0 mode required?
          */
         IterativeSchemeBase(EPMFloat a, EPMFloat b, const BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0);

         /**
          * @brief Destructor
          */
         virtual ~IterativeSchemeBase() {};

         /**
          * @brief Add a boundary condition
          *
          * @param pBC Boundary condition
          */
         virtual void addBC(SmartBC pBC) = 0;
         
      protected:
         /**
          * @brief  Iterator to current iteration
          */
         std::vector<SmartIteration>::iterator  mCurrentIt;

         /**
          * @brief \f$a\f$ coefficient of timestep operator
          */
         EPMFloat    mA;

         /**
          * @brief \f$b\f$ coefficient of timestep operator
          */
         EPMFloat    mB;

         /**
          * @brief Reference to the radial basis
          */
         const BasisType&   mrBasis;

         /**
          * @brief Vector of all the steps required for the implemented scheme
          */
         std::vector<SmartIteration>  mSteps;

         /**
          * @brief Storage for possible required intermediate values
          */
         std::vector<SmartScalarType>  mTmpScalars;

         /**
          * @brief Reset current pointer to first step
          */
         void resetIterations();

         /**
          * @brief Compute next scheme iteration
          *
          * @param rVar Input/Output variable
          * @param rNTerms New non linear terms
          */
         void doIteration(ScalarType& rVar, ScalarType& rNTerms);

      private:
   };

}

#endif // ITERATIVESCHEMEBASE_HPP
