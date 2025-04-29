/** \file ETD2RKMethod.hpp
 *  \brief Implementation of the ETD2RK method (without influence matrix)
 */

#ifndef ETD2RKMETHOD_HPP
#define ETD2RKMETHOD_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "Timestepping/TimestepParameters.hpp"
#include "Timestepping/IterativeSchemeBase.hpp"
#include "Timestepping/ETD/Operators/ETD2Operators.hpp"
#include "Timestepping/ETD/Iterations/ETD1Iteration.hpp"
#include "Timestepping/ETD/Iterations/ETD2RKTimestep.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the theta method (without influence matrix)
    */
   class ETD2RKMethod: public IterativeSchemeBase
   {
      public:
         /// Typedef from Simulation trait to local radial basis type
         typedef RadialBasisType    BasisType;

         /// Typedef from Simulation trait to local scalar type
         typedef SimulationConfig::NumericalScheme::ScalarType    ScalarType;

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
         ETD2RKMethod(EPMFloat a, EPMFloat b, const BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0);

         /**
          * @brief Destructor
          */
         virtual ~ETD2RKMethod() {};

         /**
          * @brief Add a boundary condition
          *
          * @param pBC Boundary condition
          */
         virtual void addBC(SmartBC pBC);

         /**
          * @brief Initialise the ETD2RK method
          */
         void init();
         
      protected:

         /**
          * @brief The set of ETD2 operators
          */
         ETD2Operators  mETD2;

         /**
          * @brief Update the timestep matrices after a timestep change
          */
         void updateTimeMatrices();

         /**
          * @brief Init the operators
          */
         void initOperators();

      private:

         /**
          * @brief Initialise the storage for the method
          */
         void initStorage();
   };

}

#endif // ETD2RKMETHOD_HPP
