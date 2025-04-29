/** \file NavierStokesBase.hpp
 *  \brief Base for all the Navier-Stokes equation implementations
 */

#ifndef NAVIERSTOKESBASE_HPP
#define NAVIERSTOKESBASE_HPP

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
#include "Equations/TorPolDiffusionEquation.hpp"

namespace EPMDynamo {

   /**
    * @brief General representation of the Navier-Stokes equation
    *
    * \tparam TSimTraits Traits of the simulation implementation
    * \tparam TInfluenceTraits Traits for the influence matrix step
    */
   template <typename TSimTraits, template <typename> class TInfluenceTraits> class NavierStokesBase : public TorPolDiffusionEquation<typename TSimTraits::VelType, TInfluenceTraits>  {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef SimulationConfig::TransformType    TransformType;

         /// Typedef from Simulation trait to local truncation type
         typedef SimulationConfig::NumericalScheme::ScalarType    ScalarType;

         /// Typedef for the EquationParameters type
         typedef SimulationConfig::EquationParametersType EquationParametersType;

         /**
          * @brief Constructor
          *
          * @param rV Velocity field (stored as reference)
          * @param transform Transform object (stored as reference)
          * \param tsteps Timestep parameters
          * @param params Simulation equation paramters
          */
         NavierStokesBase(typename TSimTraits::VelType &rV, TransformType &transform, TimestepParameters &tsteps, EquationParametersType &params);

         /**
          * @brief Simple empty destructor
          */
         virtual ~NavierStokesBase() {};

         /**
          * @brief Update RHS of the equation
          */
         void updateRHS() = 0;

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

   template <typename TSimTraits, template <typename> class TInfluenceTraits> inline int NavierStokesBase<TSimTraits, TInfluenceTraits>::nSSHBPacks() const
   {
      #ifdef EPMDYNAMO_RADIAL_GROUPEDCOMM
         return 3*(TSimTraits::NeedVelocity + TSimTraits::NeedVelocityCurl);
      #else
         return 3*(TSimTraits::NeedVelocity || TSimTraits::NeedVelocityCurl);
      #endif // EPMDYNAMO_RADIAL_GROUPEDCOMM
   }

   template <typename TSimTraits, template <typename> class TInfluenceTraits> inline int NavierStokesBase<TSimTraits, TInfluenceTraits>::nSHBPacks() const
   {
      #ifdef EPMDYNAMO_SH_GROUPEDCOMM
         return 3*(TSimTraits::NeedVelocity + TSimTraits::NeedVelocityCurl);
      #else
         return (TSimTraits::NeedVelocity || TSimTraits::NeedVelocityCurl);
      #endif // EPMDYNAMO_SH_GROUPEDCOMM
   }

   template <typename TSimTraits, template <typename> class TInfluenceTraits> NavierStokesBase<TSimTraits, TInfluenceTraits>::NavierStokesBase(typename TSimTraits::VelType &rV, typename NavierStokesBase<TSimTraits, TInfluenceTraits>::TransformType &transform, TimestepParameters &tsteps,  typename NavierStokesBase<TSimTraits, TInfluenceTraits>::EquationParametersType &params)
      : TorPolDiffusionEquation<typename TSimTraits::VelType, TInfluenceTraits>(rV, transform, tsteps, 1, 2, params.nsDt(), params.nsDiffusion()), mrParams(params)
   {
   }

   template <typename TSimTraits, template <typename> class TInfluenceTraits> void NavierStokesBase<TSimTraits, TInfluenceTraits>::init()
   {
      // Check that the right number of BCs have been provided
      if(this->hasAllBCs())
      {
         // Init the timesteppers (they depend on the boundary conditions)
         this->initTSteppers();
      }
      else
      {
         throw EPMException("NavierStokesBase::init", "Tried to initialise with wrong number of BCs");
      }
   }
}

#endif // NAVIERSTOKESBASE_HPP
