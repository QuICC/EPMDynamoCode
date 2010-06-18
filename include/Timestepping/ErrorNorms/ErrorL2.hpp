/** \file ErrorL2.hpp
 *  \brief Implementation of the L2 error norm computation
 */

#ifndef ERRORL2_HPP
#define ERRORL2_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the L2 error norm computation
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ErrorL2
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::ScalarType    ScalarType;

         /**
          * @brief Compute the L2 error norm
          *
          * @param rVar Input variable
          */
         static EPMFloat computeNorm(const ScalarType& rVar);

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
         ErrorL2();

         /**
          * @brief Simple empty destructor
          */
         virtual ~ErrorL2() {};
   };

   template <typename TSimType> inline EPMFloat ErrorL2<TSimType>::computeNorm(const typename ErrorL2<TSimType>::ScalarType& rVar)
   {
      return rVar.l2Norm();
   }

   template <typename TSimType> inline EPMFloat ErrorL2<TSimType>::updateNorm(const EPMFloat newNorm, const EPMFloat oldNorm)
   {
      return std::max(newNorm, oldNorm);
   }

   template <typename TSimType> ErrorL2<TSimType>::ErrorL2()
   {
   }

}

#endif // ERRORL2_HPP
