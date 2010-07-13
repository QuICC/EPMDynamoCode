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
#include "Simulations/Traits/SimulationTraits.hpp"
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
   template <typename TSimType, template <typename> class TSimTraits> class TransportBase: public ScalarDiffusionEquation<TSimType, typename TSimTraits<TSimType>::CodType >
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /**
         * @brief Constructor
         *
         * @param rC Codensity scalar (stored as reference)
         * @param transform Transform object (stored as reference)
         * \param tsteps Timestep parameters
         * @param params Simulation equation parameters
         */
         TransportBase(typename TSimTraits<TSimType>::CodType &rC, TransformType &transform, TimestepParameters &tsteps, const EquationParameters &params);

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

         /**
          * @brief Get the maximum number of backward SSH transform packs required
          *
          * The actual number of required packs depends on exact definition of equation
          */
         virtual int nSSHBPacks() const;

         /**
          * @brief Get the maximum number of backward SH transform packs required
          *
          * The actual number of required packs depends on exact definition of equation
          */
         virtual int nSHBPacks() const;
         
      protected:

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits> inline int TransportBase<TSimType, TSimTraits>::nSSHBPacks() const
   {
      #ifdef EPMDYNAMO_RADIAL_GROUPEDCOMM
         return TSimTraits<TSimType>::NeedCodensity + 2*TSimTraits<TSimType>::NeedCodensityGrad;
      #else
         return std::max(static_cast<int>(TSimTraits<TSimType>::NeedCodensity), static_cast<int>(2*TSimTraits<TSimType>::NeedCodensityGrad));
      #endif // EPMDYNAMO_RADIAL_GROUPEDCOMM
   }

   template <typename TSimType, template <typename> class TSimTraits> inline int TransportBase<TSimType, TSimTraits>::nSHBPacks() const
   {
      #ifdef EPMDYNAMO_SH_GROUPEDCOMM
         return TSimTraits<TSimType>::NeedCodensity + 3*TSimTraits<TSimType>::NeedCodensityGrad;
      #else
         return std::max(static_cast<int>(TSimTraits<TSimType>::NeedCodensity), static_cast<int>(2*TSimTraits<TSimType>::NeedCodensityGrad));
      #endif // EPMDYNAMO_SH_GROUPEDCOMM
   }

   template <typename TSimType, template <typename> class TSimTraits> TransportBase<TSimType, TSimTraits>::TransportBase(typename TSimTraits<TSimType>::CodType &rC, typename TransportBase<TSimType, TSimTraits>::TransformType &transform, TimestepParameters &tsteps, const EquationParameters &params)
      : ScalarDiffusionEquation<TSimType, typename TSimTraits<TSimType>::CodType >(rC, transform, tsteps, 1, 1.0, params.q())
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
