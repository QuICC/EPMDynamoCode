/** \file ErrorControl.hpp
 *  \brief Implementation of a specialised error control
 */

#ifndef ERRORCONTROL_HPP
#define ERRORCONTROL_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "Timestepping/ErrorNorms/ErrorL2.hpp"
#include "Timestepping/ErrorNorms/ErrorL2Summed.hpp"
#include "Timestepping/ErrorNorms/ErrorL2Max.hpp"

namespace EPMDynamo {

   /**
    * @brief This class is reponsible for providing the error norm computation routines
    *        How the norm is actually computed is defined in the timestep scheme traits
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ErrorControl
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::ScalarType    ScalarType;

         /// Typedef from Simulation trait to local truncation type
         typedef typename SimulationTraits<TSimType>::TimestepTraits::MethodTraits::ErrorNormType    ErrorNormType;

         /**
          * @brief Compute the error norm
          *
          * @param rVar Input variable
          * @param oldError Previous error norm
          */
         static EPMFloat errorNorm(const ScalarType& rVar, const EPMFloat oldError);

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

   template <typename TSimType> EPMFloat ErrorControl<TSimType>::errorNorm(const typename ErrorControl<TSimType>::ScalarType& rVar, const EPMFloat oldError)
   {
      // Create temporary storage
      EPMFloat norm = 0.0;
      EPMFloat error = 0.0;

      // Compute error norm
      norm = ErrorNormType::computeNorm(rVar);

      // Update error norm
      error = ErrorNormType::updateNorm(norm, oldError);

      return error;
   }

}

#endif // ERRORCONTROL_HPP
