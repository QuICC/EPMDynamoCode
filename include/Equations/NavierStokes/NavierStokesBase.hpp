/** \file NavierStokesBase.hpp
 *  \brief Base for all the Navier-Stokes equation implementations
 */

#ifndef NAVIERSTOKESBASE_HPP
#define NAVIERSTOKESBASE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "General/EPMTypedefs.hpp"
#include "Equations/TorPolDiffusionEquation.hpp"
#include "Equations/Parameters/EquationParameters.hpp"

namespace EPMDynamo {

   /**
    * @brief General representation of the Navier-Stokes equation
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the simulation implementation
    * \tparam TInfluenceTraits Traits for the influence matrix step
    */
   template <typename TSimType, template <typename> class TSimTraits, template <typename> class TInfluenceTraits> class NavierStokesBase : public TorPolDiffusionEquation<TSimType, typename TSimTraits<TSimType>::VelType, TInfluenceTraits>  {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::ScalarType    ScalarType;

         /**
          * @brief Constructor
          *
          * @param rV Velocity field (stored as reference)
          * @param transform Transform object (stored as reference)
          * \param tsteps Timestep parameters
          * @param params Simulation equation paramters
          */
         NavierStokesBase(typename TSimTraits<TSimType>::VelType &rV, TransformType &transform, TimestepParameters &tsteps, EquationParameters &params);

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
         EquationParameters mrParams;

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits, template <typename> class TInfluenceTraits> inline int NavierStokesBase<TSimType, TSimTraits, TInfluenceTraits>::nSSHBPacks() const
   {
      #ifdef EPMDYNAMO_SSH_GROUPEDCOMM
         return 3*(TSimTraits<TSimType>::NeedVelocity + TSimTraits<TSimType>::NeedVelocityCurl);
      #else
         return 3*(TSimTraits<TSimType>::NeedVelocity && TSimTraits<TSimType>::NeedVelocityCurl);
      #endif // EPMDYNAMO_SSH_GROUPEDCOMM
   }

   template <typename TSimType, template <typename> class TSimTraits, template <typename> class TInfluenceTraits> inline int NavierStokesBase<TSimType, TSimTraits, TInfluenceTraits>::nSHBPacks() const
   {
      #ifdef EPMDYNAMO_SH_GROUPEDCOMM
         return 3*(TSimTraits<TSimType>::NeedVelocity + TSimTraits<TSimType>::NeedVelocityCurl);
      #else
         return (TSimTraits<TSimType>::NeedVelocity && TSimTraits<TSimType>::NeedVelocityCurl);
      #endif // EPMDYNAMO_SH_GROUPEDCOMM
   }

   template <typename TSimType, template <typename> class TSimTraits, template <typename> class TInfluenceTraits> NavierStokesBase<TSimType, TSimTraits, TInfluenceTraits>::NavierStokesBase(typename TSimTraits<TSimType>::VelType &rV, typename NavierStokesBase<TSimType, TSimTraits, TInfluenceTraits>::TransformType &transform, TimestepParameters &tsteps, EquationParameters &params)
      : TorPolDiffusionEquation<TSimType, typename TSimTraits<TSimType>::VelType, TInfluenceTraits>(rV, transform, tsteps, 1, 2, params.Ro(), params.E()), mrParams(params)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits, template <typename> class TInfluenceTraits> void NavierStokesBase<TSimType, TSimTraits, TInfluenceTraits>::init()
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
