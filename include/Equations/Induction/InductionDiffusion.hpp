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
#include "Simulations/Traits/SimulationTraits.hpp"
#include "General/EPMTypedefs.hpp"
#include "Equations/Induction/InductionBase.hpp"

namespace EPMDynamo {

   /**
    * @brief General representation of the Induction equation with only diffusion included
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class InductionDiffusion : public InductionBase<TSimType, TSimTraits>
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
         InductionDiffusion(typename TSimTraits<TSimType>::MagType &rB, TransformType &transform, TimestepParameters &tsteps);

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

   template <typename TSimType, template <typename> class TSimTraits> InductionDiffusion<TSimType, TSimTraits>::InductionDiffusion(typename TSimTraits<TSimType>::MagType &rB, typename InductionDiffusion<TSimType, TSimTraits>::TransformType &transform, TimestepParameters &tsteps)
      : InductionBase<TSimType, TSimTraits>(rB, transform, tsteps)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void InductionDiffusion<TSimType, TSimTraits>::updateRTP(const int step)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void InductionDiffusion<TSimType, TSimTraits>::updateRHS()
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void InductionDiffusion<TSimType, TSimTraits>::transformRHS(const int step)
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
