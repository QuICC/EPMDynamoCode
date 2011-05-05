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
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class InductionBase : public TorPolDiffusionEquation<TSimType, typename TSimTraits<TSimType>::MagType>
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
         InductionBase(typename TSimTraits<TSimType>::MagType &rB, TransformType &transform, TimestepParameters &tsteps, EquationParametersType &params);

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

   template <typename TSimType, template <typename> class TSimTraits> inline int InductionBase<TSimType, TSimTraits>::nSSHBPacks() const
   {
      #ifdef EPMDYNAMO_RADIAL_GROUPEDCOMM
         return 3*(TSimTraits<TSimType>::NeedMagnetic + TSimTraits<TSimType>::NeedMagneticCurl);
      #else
         return 3*(TSimTraits<TSimType>::NeedMagnetic || TSimTraits<TSimType>::NeedMagneticCurl);
      #endif // EPMDYNAMO_RADIAL_GROUPEDCOMM
   }

   template <typename TSimType, template <typename> class TSimTraits> inline int InductionBase<TSimType, TSimTraits>::nSHBPacks() const
   {
      #ifdef EPMDYNAMO_SH_GROUPEDCOMM
         return 3*(TSimTraits<TSimType>::NeedMagnetic + TSimTraits<TSimType>::NeedMagneticCurl);
      #else
         return (TSimTraits<TSimType>::NeedMagnetic || TSimTraits<TSimType>::NeedMagneticCurl);
      #endif // EPMDYNAMO_SH_GROUPEDCOMM
   }

   template <typename TSimType, template <typename> class TSimTraits> InductionBase<TSimType, TSimTraits>::InductionBase(typename TSimTraits<TSimType>::MagType &rB, typename InductionBase<TSimType, TSimTraits>::TransformType &transform, TimestepParameters &tsteps, typename InductionBase<TSimType, TSimTraits>::EquationParametersType &params)
      : TorPolDiffusionEquation<TSimType, typename TSimTraits<TSimType>::MagType>(rB, transform, tsteps, 1, 1, params.indDt(), params.indDiffusion()), mrParams(params)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void InductionBase<TSimType, TSimTraits>::init()
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
