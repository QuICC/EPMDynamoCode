/** \file ErrorModeL2Max.hpp
 *  \brief Implementation of the relative L2 max error norm per harmonic mode computation
 */

#ifndef ERRORMODEL2MAX_HPP
#define ERRORMODEL2MAX_HPP

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
   class ErrorModeL2Max
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
         ErrorModeL2Max();

         /**
          * @brief Simple empty destructor
          */
         virtual ~ErrorModeL2Max() {}; 
   };

   inline EPMFloat ErrorModeL2Max::updateNorm(const EPMFloat newNorm, const EPMFloat oldNorm)
   {
      return std::max(newNorm, oldNorm);
   }

}

#endif // ERRORMODEL2MAX_HPP
