/** \file InductionBase.hpp
 *  \brief Base for all the possible induction equations
 */

#ifndef INDUCTIONBASE_HPP
#define INDUCTIONBASE_HPP

// System includes
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
    * @brief General representation of the Induction equation
    *
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimTraits> class InductionBase : public TorPolDiffusionEquation<typename TSimTraits::MagType>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef SimulationConfig::TransformType    TransformType;

         /// Typedef for the EquationParameters type
         typedef SimulationConfig::EquationParametersType EquationParametersType;

         /**
          * @brief Constructor
          *
          * @param rB Magnetic field (stored as reference)
          * @param transform Transform object (stored as reference)
          * @param tsteps Timestep parameters
          * @param params Equation parameters
          */
         InductionBase(typename TSimTraits::MagType &rB, TransformType &transform, TimestepParameters &tsteps, EquationParametersType &params);

         /**
          * @brief Simple empty destructor
          */
         virtual ~InductionBase() {};

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

   template <typename TSimTraits> inline int InductionBase<TSimTraits>::nSSHBPacks() const
   {
      #ifdef EPMDYNAMO_RADIAL_GROUPEDCOMM
         return 3*(TSimTraits::NeedMagnetic + TSimTraits::NeedMagneticCurl);
      #else
         return 3*(TSimTraits::NeedMagnetic || TSimTraits::NeedMagneticCurl);
      #endif // EPMDYNAMO_RADIAL_GROUPEDCOMM
   }

   template <typename TSimTraits> inline int InductionBase<TSimTraits>::nSHBPacks() const
   {
      #ifdef EPMDYNAMO_SH_GROUPEDCOMM
         return 3*(TSimTraits::NeedMagnetic + TSimTraits::NeedMagneticCurl);
      #else
         return (TSimTraits::NeedMagnetic || TSimTraits::NeedMagneticCurl);
      #endif // EPMDYNAMO_SH_GROUPEDCOMM
   }

   template <typename TSimTraits> InductionBase<TSimTraits>::InductionBase(typename TSimTraits::MagType &rB, typename InductionBase<TSimTraits>::TransformType &transform, TimestepParameters &tsteps, typename InductionBase<TSimTraits>::EquationParametersType &params)
      : TorPolDiffusionEquation<typename TSimTraits::MagType>(rB, transform, tsteps, 1, 1, params.indDt(), params.indDiffusion()), mrParams(params)
   {
   }

   template <typename TSimTraits> void InductionBase<TSimTraits>::init()
   {
      // Check that the right number of BCs have been provided
      if(this->hasAllBCs())
      {
         // Init the timesteppers (they depend on the boundary conditions)
         this->initTSteppers();
      }
      else
      {
         throw EPMException("InductionBase::init", "Tried to initialise with wrong number of BCs");
      }
   }
}

#endif // INDUCTIONBASE_HPP
