/** \file ErrorCoeffL2Summed.hpp
 *  \brief Implementation of the relative L2 max error norm per radial coefficient computation
 */

#ifndef ERRORCOEFFL2SUMMED_HPP
#define ERRORCOEFFL2SUMMED_HPP

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
   class ErrorCoeffL2Summed
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
         ErrorCoeffL2Summed();

         /**
          * @brief Simple empty destructor
          */
         virtual ~ErrorCoeffL2Summed() {}; 
   };

   inline EPMFloat ErrorCoeffL2Summed::updateNorm(const EPMFloat newNorm, const EPMFloat oldNorm)
   {
      return std::max(newNorm, oldNorm);
   }

}

#endif // ERRORCOEFFL2SUMMED_HPP
