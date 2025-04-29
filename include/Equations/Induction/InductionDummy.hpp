/** \file InductionDummy.hpp
 *  \brief Implementatio of the DUMMY induction equation
 */

#ifndef INDUCTIONDUMMY_HPP
#define INDUCTIONDUMMY_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Config/SimulationConfig.hpp"
#include "General/EPMTypedefs.hpp"
#include "Equations/Induction/InductionBase.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of the full DUMMY Induction equation
    *
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimTraits> class InductionDummy : public InductionBase<TSimTraits>
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
          * @param rV Velocity field (stored as reference)
          * @param transform Transform object (stored as reference)
          * @param tsteps Timestep parameters
          * @param params Equation parameters
          */
         InductionDummy(typename TSimTraits::MagType &rB, typename TSimTraits::VelType &rV, TransformType &transform, TimestepParameters &tsteps, EquationParametersType &params);

         /**
          * @brief Simple empty destructor
          */
         virtual ~InductionDummy() {};

         /**
          * @brief Update RTP values of the equation
          *
          * \param step Current step in a multistep transform
          */
         void updateRTP(const int step);

         /**
          * @brief Update RHS of the equation
          */
         void updateRHS();

         /**
          * @brief Transform RHS of the equation
          *
          * \param step Current step in a multistep transform
          */
         void transformRHS(const int step);

         /*
          * @brief Dummy timestepping to trigger backward transforms
          */
         void timestep();
         
      protected:

         /**
          * @brief Const Reference variable to the velocity field
          */
         typename TSimTraits::VelType&  mrV;

      private:
   };

   template <typename TSimTraits> InductionDummy<TSimTraits>::InductionDummy(typename TSimTraits::MagType &rB, typename TSimTraits::VelType &rV, typename InductionDummy<TSimTraits>::TransformType &transform, TimestepParameters &tsteps,  typename InductionDummy<TSimTraits>::EquationParametersType &params)
      : InductionBase<TSimTraits>(rB, transform, tsteps, params), mrV(rV)
   {
   }

   template <typename TSimTraits> void InductionDummy<TSimTraits>::updateRTP(const int step)
   {
      // Update real space values for Magnetic field
      this->mrX.rOc().transform(step);

      // Update real space values for velocity field
      this->mrV.rOc().transform(step);
   }

   template <typename TSimTraits> void InductionDummy<TSimTraits>::updateRHS()
   {
      // Compute cross product \f$\vec{u}\times\vec{B}\f$
      this->mrV.oc().rtp().template cross<0>(this->mNTerms.rOc().rRTP(), this->mrX.oc().rtp(), this->mrParams.indAdvection());
   }

   template <typename TSimTraits> void InductionDummy<TSimTraits>::transformRHS(const int step)
   {
      // Transform Non Linear terms to spectral space from mNTerms values
      this->transformNTerms(this->mNTerms.rOc().rPerturbation().rPol(), this->mNTerms.rOc().rPerturbation().rTor());
   }

   template <typename TSimTraits> void InductionDummy<TSimTraits>::timestep()
   {
      this->mrX.rOc().rPerturbation().rTor();
      this->mrX.rOc().rPerturbation().rPol();
   }
}

#endif // INDUCTIONDUMMY_HPP
