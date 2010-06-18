/** \file PCErrorControl.hpp
 *  \brief Implementation of a specialised error control for predictor/corrector schemes
 */

#ifndef PCERRORCONTROL_HPP
#define PCERRORCONTROL_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/SimulationTraits.hpp"
#include "Timestepping/ErrorNorms/ErrorL2.hpp"
#include "Timestepping/ErrorNorms/ErrorL2Summed.hpp"
#include "Timestepping/ErrorNorms/ErrorL2Max.hpp"

namespace EPMDynamo {

   /**
    * @brief This class is reponsible for providing the error norm computation routines
    *        How the norm is actually computed is defined in the timestep scheme traits
    *
    * \tparam TSimType Type of the simulation
    *
    * \bug Needs a big cleaning and rewritting!
    */
   template <typename TSimType> class PCErrorControl
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::ScalarType    ScalarType;

         /// Typedef from Simulation trait to local truncation type
         typedef typename SimulationTraits<TSimType>::TimestepTraits::ErrorNormType    ErrorNormType;

         /**
          * @brief Compute the error norm
          *
          * @param rVar Input variable
          * @param oldError Previous error norm
          */
         static DynamoFloat errorNorm(const ScalarType& rVar, const DynamoFloat oldError);

      protected:

      private:
         /**
          * @brief Constructor
          */
         PCErrorControl() {};

         /**
          * @brief Simple empty destructor
          */
         virtual ~PCErrorControl() {};
   };

   template <typename TSimType> DynamoFloat PCErrorControl<TSimType>::errorNorm(const typename PCErrorControl<TSimType>::ScalarType& rVar, const DynamoFloat oldError)
   {
      // Create temporary storage
      DynamoFloat norm = 0.0;
      DynamoFloat error = 0.0;

      // Compute error norm
      norm = ErrorNormType::computeNorm(rVar);

      // Update error norm
      error = ErrorNormType::updateNorm(norm, oldError);

      return error;
   }

}

#endif // PCERRORCONTROL_HPP
