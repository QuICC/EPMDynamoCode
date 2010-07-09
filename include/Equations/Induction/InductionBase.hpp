/** \file InductionBase.hpp
 *  \brief Base for all the possible induction equations
 */

#ifndef INDUCTIONBASE_HPP
#define INDUCTIONBASE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
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
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /**
          * @brief Constructor
          *
          * @param rB Magnetic field (stored as reference)
          * @param transform Transform object (stored as reference)
          * \param tsteps Timestep parameters
          */
         InductionBase(typename TSimTraits<TSimType>::MagType &rB, TransformType &transform, TimestepParameters &tsteps);

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

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits> inline int InductionBase<TSimType, TSimTraits>::nSSHBPacks() const
   {
      #ifdef EPMDYNAMO_SSH_GROUPEDCOMM
         return 3*(TSimTraits<TSimType>::NeedMagnetic + TSimTraits<TSimType>::NeedMagneticCurl);
      #else
         return 3*(TSimTraits<TSimType>::NeedMagnetic || TSimTraits<TSimType>::NeedMagneticCurl);
      #endif // EPMDYNAMO_SSH_GROUPEDCOMM
   }

   template <typename TSimType, template <typename> class TSimTraits> inline int InductionBase<TSimType, TSimTraits>::nSHBPacks() const
   {
      #ifdef EPMDYNAMO_SH_GROUPEDCOMM
         return 3*(TSimTraits<TSimType>::NeedMagnetic + TSimTraits<TSimType>::NeedMagneticCurl);
      #else
         return (TSimTraits<TSimType>::NeedMagnetic || TSimTraits<TSimType>::NeedMagneticCurl);
      #endif // EPMDYNAMO_SH_GROUPEDCOMM
   }

   template <typename TSimType, template <typename> class TSimTraits> InductionBase<TSimType, TSimTraits>::InductionBase(typename TSimTraits<TSimType>::MagType &rB, typename InductionBase<TSimType, TSimTraits>::TransformType &transform, TimestepParameters &tsteps)
      : TorPolDiffusionEquation<TSimType, typename TSimTraits<TSimType>::MagType>(rB, transform, tsteps, 1, 1, 1.0, 1.0)
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
