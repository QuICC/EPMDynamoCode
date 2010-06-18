/** \file ErrorL2Summed.hpp
 *  \brief Implementation of the L2 summed error norm computation
 */

#ifndef ERRORL2SUMMED_HPP
#define ERRORL2SUMMED_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the L2 summed error norm computation
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ErrorL2Summed
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::ScalarType    ScalarType;

         /**
          * @brief Compute the L2 summed error norm
          *
          * @param rVar Input variable
          */
         static EPMFloat computeNorm(const ScalarType& rVar);

         /**
          * @brief Add new to old value
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
         ErrorL2Summed();

         /**
          * @brief Simple empty destructor
          */
         virtual ~ErrorL2Summed() {};
   };

   template <typename TSimType> inline EPMFloat ErrorL2Summed<TSimType>::computeNorm(const typename ErrorL2Summed<TSimType>::ScalarType& rVar)
   {
      //return rVar.l2Norm();
      return rVar.modalL2Norm();
   }

   template <typename TSimType> inline EPMFloat ErrorL2Summed<TSimType>::updateNorm(const EPMFloat newNorm, const EPMFloat oldNorm)
   {
      return newNorm + oldNorm;
   }

   template <typename TSimType> ErrorL2Summed<TSimType>::ErrorL2Summed()
   {
   }

}

#endif // ERRORL2SUMMED_HPP
