/** \file ETD2RKInfluenceMethod.hpp
 *  \brief Implementation of the ETD2RK method (with influence matrix)
 */

#ifndef ETD2RKINFLUENCEMETHOD_HPP
#define ETD2RKINFLUENCEMETHOD_HPP

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
#include "Timestepping/ETD/ETD2RK/ETD2RKMethod.hpp"
#include "Timestepping/ETD/ETDInfluenceKernel.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the theta method (without influence matrix)
    */
   class ETD2RKInfluenceMethod: public ETD2RKMethod
   {
      public:
         /// Typedef from Simulation trait to local radial basis type
         typedef SimulationConfig::NumericalScheme::RadialBasisType    BasisType;

         /// Typedef from Simulation trait to local scalar type
         typedef SimulationConfig::NumericalScheme::ScalarType    ScalarType;

         /// Typedef for a smart ETD influence kernel
         typedef EPMSHARED_PTR<ETDInfluenceKernel > SmartInfluenceKernel;

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
         ETD2RKInfluenceMethod(EPMFloat a, EPMFloat b, const BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0);

         /**
          * @brief Destructor
          */
         virtual ~ETD2RKInfluenceMethod() {};

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
          * @brief Initialise the ETD2 influence matrix method
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
         InfluenceMatrix mInfluence;

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

}

#endif // ETD2RKINFLUENCEMETHOD_HPP
