/** \file ETD2RKTimestep.hpp
 *  \brief Implementation of the ETD2RK timestep step
 */

#ifndef ETD2RKTIMESTEP_HPP
#define ETD2RKTIMESTEP_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "Timestepping/SchemeIteration.hpp"
#include "Timestepping/ETD/ETDSchemeTraits.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of the ETD2RK timestep step
    *
    * This class performs the actual timestep comuputation after the intermediate
    * computations are performed.
    */
   class ETD2RKTimestep: public SchemeIteration
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef SimulationConfig::NumericalScheme::ScalarType    ScalarType;

         /// typedef for a pointer to a scalar type
         typedef EPMSHARED_PTR<ScalarType> SmartScalarType;

         /// Typedef for a smart pointer to ETDOperators
         typedef EPMSHARED_PTR<ETDSchemeTraits::Operators> SmartETDOperators;

         /**
          * @brief Constructor
          *
          * @param nFactor Non linear terms multiplicative factor
          * @param pOldN Pointer to the old non linear terms
          * @param pOpM2 Pointer to the \f$M_2\f$ operator
          */
         ETD2RKTimestep(EPMFloat nFactor, SmartScalarType pOldN, SmartETDOperators pOpM2);

         /**
          * @brief Destructor
          */
         virtual ~ETD2RKTimestep() {};

         /**
          * @brief Prepare the computation of the implemented iteration
          *
          * @param rVar Input/Output variable
          * @param rNTerms New non linear terms
          */
         virtual void prepare(ScalarType &rVar, ScalarType &rNTerms);

         /**
          * @brief Compute the implemented iteration
          *
          * @param rVar Input/Output variable
          * @param rNTerms New non linear terms
          */
         virtual void compute(ScalarType &rVar, ScalarType &rNTerms);
         
      protected:
         /**
          * @brief multiplicative factor required for the non linear terms
          */
         EPMFloat mNFactor;

         /**
          * @brief Storage for the old NTerms
          */
         SmartScalarType   mpOldN;

         /**
          * @brief ETD operator M2
          */
         SmartETDOperators  mpOpM2;

      private:
   };

}

#endif // ETD2RKTIMESTEP_HPP
