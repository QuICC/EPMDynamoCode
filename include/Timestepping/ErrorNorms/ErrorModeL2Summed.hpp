/** \file ErrorModeL2Summed.hpp
 *  \brief Implementation of the relative L2 max error norm per harmonic mode computation
 */

#ifndef ERRORMODEL2SUMMED_HPP
#define ERRORMODEL2SUMMED_HPP

// Configuration includes
//
#include "Config/Parallelisation.h"
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * \brief Implementation of the L2 max error norm computation
    */
   class ErrorModeL2Summed
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef SimulationConfig::NumericalScheme::ScalarType    ScalarType;

         /**
          * @brief Compute the Max error norm
          *
          * @param rVar Input variable
          * @param rRef Reference variable for relative error
          */
         static EPMFloat computeNorm(const ScalarType& rVar, const ScalarType& rRef);

         /**
          * @brief Compare new and old values and keep maximum
          *
          * @param newNorm New error norm value
          * @param oldNorm Old error norm value
          */
         static EPMFloat updateNorm(const EPMFloat newNorm, const EPMFloat oldNorm);

      protected:

      private:
         /**
          * @brief Constructor
          */
         ErrorModeL2Summed();

         /**
          * @brief Simple empty destructor
          */
         virtual ~ErrorModeL2Summed() {}; 
   };

   inline EPMFloat ErrorModeL2Summed::updateNorm(const EPMFloat newNorm, const EPMFloat oldNorm)
   {
      return std::max(newNorm, oldNorm);
   }

}

#endif // ERRORMODEL2SUMMED_HPP
