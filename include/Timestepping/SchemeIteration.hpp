/** \file SchemeIteration.hpp
 *  \brief General implementation of an iterative scheme iteration (step)
 */

#ifndef SCHEMEITERATION_HPP
#define SCHEMEITERATION_HPP

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

namespace EPMDynamo {

   /**
    * \brief General implementation of an iterative scheme iteration (step)
    */
   class SchemeIteration
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef SimulationConfig::NumericalScheme::ScalarType    ScalarType;

         /**
          * @brief Constructor
          */
         SchemeIteration();

         /**
          * @brief Destructor
          */
         virtual ~SchemeIteration() {};

         /**
          * @brief Prepare the computation of the iteration
          *
          * @param rVar Input/Output variable
          * @param rNTerms New non linear terms
          */
         virtual void prepare(ScalarType &rVar, ScalarType &rNTerms){};

         /**
          * @brief Compute the implemented iteration
          *
          * @param rVar Input/Output variable
          * @param rNTerms New non linear terms
          */
         virtual void compute(ScalarType &rVar, ScalarType &rNTerms) = 0;

         /**
          * @brief Does this iteration provide timestep error information
          */
         bool providesError() const;

         /**
          * @brief Does this iteration provide timestep error information
          */
         bool doNextIteration() const;
         
      protected:

      private:
   };

   inline bool SchemeIteration::providesError() const
   {
      return false;
   }

   inline bool SchemeIteration::doNextIteration() const
   {
      return true;
   }
}

#endif // SCHEMEITERATION_HPP
