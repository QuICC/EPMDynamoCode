/** \file TransportBase.hpp
 *  \brief Base for all the transport equations
 */

#ifndef TRANSPORTBASE_HPP
#define TRANSPORTBASE_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "General/EPMException.hpp"
#include "Equations/ScalarDiffusionEquation.hpp"

namespace EPMDynamo {

   /**
    * @brief General representation of the Transport equation
    *
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimTraits> class TransportBase: public ScalarDiffusionEquation<typename TSimTraits::CodType >
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef SimulationConfig::TransformType    TransformType;

         /// Typedef for the EquationParameters type
         typedef SimulationConfig::EquationParametersType EquationParametersType;

         /**
         * @brief Constructor
         *
         * @param rC Codensity scalar (stored as reference)
         * @param transform Transform object (stored as reference)
         * \param tsteps Timestep parameters
         * @param params Simulation equation parameters
         */
         TransportBase(typename TSimTraits::CodType &rC, TransformType &transform, TimestepParameters &tsteps, const EquationParametersType &params);

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

         /**
          * @brief Reference to the equation parameters object
          */
         EquationParametersType mrParams;

      private:
   };

   template <typename TSimTraits> inline int TransportBase<TSimTraits>::nSSHBPacks() const
   {
      #ifdef EPMDYNAMO_RADIAL_GROUPEDCOMM
         return TSimTraits::NeedCodensity + 2*TSimTraits::NeedCodensityGrad;
      #else
         return std::max(static_cast<int>(TSimTraits::NeedCodensity), static_cast<int>(2*TSimTraits::NeedCodensityGrad));
      #endif // EPMDYNAMO_RADIAL_GROUPEDCOMM
   }

   template <typename TSimTraits> inline int TransportBase<TSimTraits>::nSHBPacks() const
   {
      #ifdef EPMDYNAMO_SH_GROUPEDCOMM
         return TSimTraits::NeedCodensity + 3*TSimTraits::NeedCodensityGrad;
      #else
         return std::max(static_cast<int>(TSimTraits::NeedCodensity), static_cast<int>(2*TSimTraits::NeedCodensityGrad));
      #endif // EPMDYNAMO_SH_GROUPEDCOMM
   }

   template <typename TSimTraits> TransportBase<TSimTraits>::TransportBase(typename TSimTraits::CodType &rC, typename TransportBase<TSimTraits>::TransformType &transform, TimestepParameters &tsteps, const  typename TransportBase<TSimTraits>::EquationParametersType &params)
      : ScalarDiffusionEquation<typename TSimTraits::CodType >(rC, transform, tsteps, 1, params.tptDt(), params.tptDiffusion()), mrParams(params)
   {
   }

   template <typename TSimTraits> void TransportBase<TSimTraits>::init()
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
