/** \file ErrorControl.hpp
 *  \brief Implementation of a specialised error control
 */

#ifndef ERRORCONTROL_HPP
#define ERRORCONTROL_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"
#include "Config/NumericalSchemeInc.hpp"
#include "Config/TimesteppingInc.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "Timestepping/ErrorNorms/ErrorL2.hpp"
#include "Timestepping/ErrorNorms/ErrorCoeffL2Summed.hpp"
#include "Timestepping/ErrorNorms/ErrorModeL2Summed.hpp"
#include "Timestepping/ErrorNorms/ErrorCoeffL2Max.hpp"
#include "Timestepping/ErrorNorms/ErrorModeL2Max.hpp"

namespace EPMDynamo {

   /**
    * @brief This class is reponsible for providing the error norm computation routines
    *        How the norm is actually computed is defined in the timestep scheme traits
    */
   class ErrorControl
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef SimulationConfig::NumericalScheme::ScalarType    ScalarType;

         /// Typedef from Simulation trait to local truncation type
         typedef SimulationConfig::TimestepTraits::MethodTraits::ErrorNormType    ErrorNormType;

         /**
          * @brief Compute the error norm
          *
          * @param rVar Input variable
          * @param rRef Reference variable for relative error
          * @param oldError Previous error norm
          */
         static EPMFloat errorNorm(const ScalarType& rVar, const ScalarType& rRef, const EPMFloat oldError);

      protected:

      private:
         /**
          * @brief Constructor
          */
         ErrorControl() {};

         /**
          * @brief Simple empty destructor
          */
         virtual ~ErrorControl() {};
   };

}

#endif // ERRORCONTROL_HPP
