/** \file TransportBase.hpp
 *  \brief Base for all the transport equations
 */

#ifndef TRANSPORTBASE_HPP
#define TRANSPORTBASE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/SimulationTraits.hpp"
#include "General/EPMTypedefs.hpp"
#include "General/EPMException.hpp"
#include "Equations/Parameters/EquationParameters.hpp"
#include "Equations/ScalarDiffusionEquation.hpp"

namespace EPMDynamo {

   /**
    * @brief General representation of the Transport equation
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class TransportBase: public ScalarDiffusionEquation<TSimType, typename TSimTraits<TSim>::CodType >
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSim>::TransformType    TransformType;

         /**
         * @brief Constructor
         *
         * @param rC Codensity scalar (stored as reference)
         * @param transform Transform object (stored as reference)
         * \param tsteps Timestep parameters
         * @param params Simulation equation parameters
         */
         TransportBase(typename TSimTraits<TSim>::CodType &rC, TransformType &transform, TimestepParameters &tsteps, const EquationParameters &params);

         /**
         * @brief Simple empty destructor
         */
         virtual ~TransportBase() {};

         /**
          * @brief Update RHS of the equation
          */
         virtual void updateRHS() = 0;

         /**
          * @brief Initialise the equation
          */
         void init();
         
      protected:

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits> TransportBase<TSimType, TSimTraits>::TransportBase(typename TSimTraits<TSim>::CodType &rC, typename TransportBase<TSimType, TSimTraits>::TransformType &transform, TimestepParameters &tsteps, const EquationParameters &params)
      : ScalarDiffusionEquation<TSimType, typename TSimTraits<TSim>::CodType >(rC, transform, tsteps, 1, 1.0, params.q())
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void TransportBase<TSimType, TSimTraits>::init()
   {
      // Check that the right number of BCs have been provided
      if(this->hasAllBCs())
      {
         // Init the timesteppers (they depend on the boundary conditions)
         this->initTStepper();
      }
      else
      {
         throw EPMException("TransportBase::init", "Tried to initialise with wrong number of BCs");
      }
   }

}

#endif // TRANSPORTBASE_HPP
