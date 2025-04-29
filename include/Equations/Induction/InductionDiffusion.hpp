/** \file InductionDiffusion.hpp
 *  \brief Implementatio of the magnetic diffusion equation
 */

#ifndef INDUCTIONDIFFUSION_HPP
#define INDUCTIONDIFFUSION_HPP

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
    * @brief General representation of the Induction equation with only diffusion included
    *
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimTraits> class InductionDiffusion : public InductionBase<TSimTraits>
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
         InductionDiffusion(typename TSimTraits::MagType &rB, TransformType &transform, TimestepParameters &tsteps, EquationParametersType &params);

         /**
          * @brief Simple empty destructor
          */
         virtual ~InductionDiffusion() {};

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
         
      protected:

      private:
   };

   template <typename TSimTraits> InductionDiffusion<TSimTraits>::InductionDiffusion(typename TSimTraits::MagType &rB, typename InductionDiffusion<TSimTraits>::TransformType &transform, TimestepParameters &tsteps, typename InductionDiffusion<TSimTraits>::EquationParametersType &params)
      : InductionBase<TSimTraits>(rB, transform, tsteps, params)
   {
   }

   template <typename TSimTraits> void InductionDiffusion<TSimTraits>::updateRTP(const int step)
   {
   }

   template <typename TSimTraits> void InductionDiffusion<TSimTraits>::updateRHS()
   {
   }

   template <typename TSimTraits> void InductionDiffusion<TSimTraits>::transformRHS(const int step)
   {
      if(step == 0)
      {
         // Set the non linear terms to zero
         int nL = this->mNTerms.oc().perturbation().nL();
         int l0 = this->mNTerms.oc().perturbation().pol().minL();

         for(int l = l0; l < nL; ++l)
         {
            this->mNTerms.rOc().rPerturbation().rTor().rLShell(l).setConstant(0.0);
            this->mNTerms.rOc().rPerturbation().rPol().rLShell(l).setConstant(0.0);
         }
      }
   }
}

#endif // INDUCTIONDIFFUSION_HPP
