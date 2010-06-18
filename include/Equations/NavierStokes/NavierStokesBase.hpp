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
    * \tparam TTStepPol NEED TO CHANGE THIS
    *
    * \bug Need to change the implementation of the poisson solver
    * \bug Poisson solver has not yet been moved into new tree
    */
   template <typename TSimType, template <typename> class TSimTraits, template <typename> class TTStepPol> class NavierStokesBase : public TorPolDiffusionEquation<TSimType, typename TSimTraits<TSimType>::VelType, TTStepPol>  {
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
         
      protected:

         /**
          * @brief Reference to the equation parameters object
          */
         EquationParameters mrParams;

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits, template <typename> class TTStepPol> NavierStokesBase<TSimType, TSimTraits, TTStepPol>::NavierStokesBase(typename TSimTraits<TSimType>::VelType &rV, typename NavierStokesBase<TSimType, TSimTraits, TTStepPol>::TransformType &transform, TimestepParameters &tsteps, EquationParameters &params)
      : TorPolDiffusionEquation<TSimType, typename TSimTraits<TSimType>::VelType, TTStepPol>(rV, transform, tsteps, 1, 1, params.Ro(), params.E()), mrParams(params)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits, template <typename> class TTStepPol> void NavierStokesBase<TSimType, TSimTraits, TTStepPol>::init()
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
