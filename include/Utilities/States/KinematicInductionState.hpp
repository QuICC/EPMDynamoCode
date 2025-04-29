/** \file KinematicInductionState.hpp
 *  \brief Initial state definition(s) for kinematic induction simulation(s)
 */

#ifndef KINEMATICINDUCTIONSTATE_HPP
#define KINEMATICINDUCTIONSTATE_HPP

// Configuration includes
//
#include "Config/Parallelisation.h"

// System includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"

namespace EPMDynamo {

   /**
    * @brief Traits of kinematic induction state generator
    */
   class KinematicStateTraits
   {
      public:
         /// Requires RTP Codensity computations
         static const bool UseRTPCodensity = false;

         /// Requires spectral Codensity computations
         static const bool UseSpecCodensity = false;

         /// Requires Codensity gradient computations
         static const bool UseCodensityGrad = false;

         /// Requires RTP Magnetic computations
         static const bool UseRTPMagnetic = false;

         /// Requires spectral Magnetic computations
         static const bool UseSpecMagnetic = true;

         /// Requires Magnetic curl computations
         static const bool UseMagneticCurl = false;

         /// Requires RTP Velocity computations
         static const bool UseRTPVelocity = false;

         /// Requires Velocity computations
         static const bool UseSpecVelocity = true;

         /// Requires Velocity curl computations
         static const bool UseVelocityCurl = false;
   };

   /**
    * \brief Initial state definition(s) for the dynamo benchmark(s)
    *
    * \epmBug State generation is not aware of parallelisation
    */
   template <typename TGenTraits> class KinematicInductionState
   {
      public:
         /// Typdef for the StateTraits type
         typedef KinematicStateTraits  StateTraits;

         /**
          * @brief Velocity field case
          */
         static const int VELOCITY_CASE;

         /**
          * @brief Perturbation amplitude
          */
         static const EPMFloat PERTURBATION_AMPLITUDE;

         /**
          * @brief Ratio of L spectrum to be perturbed
          */
         static const int PERTURBATION_LRATIO;

         /**
          * @brief Ratio of radial spectrum to be perturbed
          */
         static const int PERTURBATION_NRATIO;

         /// Typedef for the codensity type
         typedef typename TGenTraits::CodType  Codensity;

         /// Typedef for the magnetic type
         typedef typename TGenTraits::MagType  Magnetic;

         /// Typedef for the velocity type
         typedef typename TGenTraits::VelType  Velocity;

         /**
          * @brief Set RTP value for the codensity scalar
          */
         static void setRTPCodensity(Codensity &codC);

         /**
          * @brief Set RTP value for the magnetic field
          */
         static void setRTPMagnetic(Magnetic &magB);

         /**
          * @brief Set RTP value for the velocity field
          */
         static void setRTPVelocity(Velocity &velV);

         /**
          * @brief Set spectral coefficients for the codensity scalar
          */
         static void setSpecCodensity(Codensity &codC);

         /**
          * @brief Set spectral coefficients for the velocity field
          */
         static void setSpecMagnetic(Magnetic &magB);

         /**
          * @brief Set spectral coefficients for the magnetic field
          */
         static void setSpecVelocity(Velocity &velV);

      private:
         /**
          * @brief Private constructor
          */
         KinematicInductionState();

         /**
          * @brief Destructor
          */
         virtual ~KinematicInductionState() {};
   };

   template <typename TGenTraits> const int KinematicInductionState<TGenTraits>::VELOCITY_CASE = 20;

   template <typename TGenTraits> const EPMFloat KinematicInductionState<TGenTraits>::PERTURBATION_AMPLITUDE = 1.0e-10;

   template <typename TGenTraits> const int KinematicInductionState<TGenTraits>::PERTURBATION_LRATIO = 2;

   template <typename TGenTraits> const int KinematicInductionState<TGenTraits>::PERTURBATION_NRATIO = 4;

   template <typename TGenTraits> void KinematicInductionState<TGenTraits>::setRTPCodensity(typename KinematicInductionState<TGenTraits>::Codensity &codC)
   {
   }

   template <typename TGenTraits> void KinematicInductionState<TGenTraits>::setRTPMagnetic(typename KinematicInductionState<TGenTraits>::Magnetic &magB)
   {
   }

   template <typename TGenTraits> void KinematicInductionState<TGenTraits>::setRTPVelocity(typename KinematicInductionState<TGenTraits>::Velocity &velV)
   {
   }

   template <typename TGenTraits> void KinematicInductionState<TGenTraits>::setSpecCodensity(typename KinematicInductionState<TGenTraits>::Codensity &codC)
   {
   }

   template <typename TGenTraits> void KinematicInductionState<TGenTraits>::setSpecMagnetic(typename KinematicInductionState<TGenTraits>::Magnetic &magB)
   {
      SmartTruncation pTrunc = magB.oc().trunc();

      for(int l=1; l < pTrunc->local()->spec()->nL()/RandomState<TGenTraits>::PERTURBATION_LRATIO; ++l)
      {
         // Set some perturbation random energy in Toroidal component
         magB.rOc().rPerturbation().rTor().rLShell(l).block(0,0,pTrunc->sim()->rad()->nN()/RandomState<TGenTraits>::PERTURBATION_NRATIO, std::max(pTrunc->local()->spec()->nM(l),1)).setRandom();
         magB.rOc().rPerturbation().rTor().rLShell(l) *= EPMComplex(RandomState<TGenTraits>::PERTURBATION_AMPLITUDE,0.0);

         // Set some perturbation random energy in Poloidal component
         magB.rOc().rPerturbation().rPol().rLShell(l).block(0,0,pTrunc->sim()->rad()->nN()/RandomState<TGenTraits>::PERTURBATION_NRATIO, std::max(pTrunc->local()->spec()->nM(l),1)).setRandom();
         magB.rOc().rPerturbation().rPol().rLShell(l) *= EPMComplex(RandomState<TGenTraits>::PERTURBATION_AMPLITUDE,0.0);

         // Make sure the m=0 imaginary part is zero!
         for(int n=0; n < pTrunc->sim()->rad()->nN(); ++n)
         {
            magB.rOc().rPerturbation().rTor().rLShell(l).col(0)(n).imag() = 0.0;
            magB.rOc().rPerturbation().rPol().rLShell(l).col(0)(n).imag() = 0.0;
         }
      }
   }

   template <typename TGenTraits> void KinematicInductionState<TGenTraits>::setSpecVelocity(typename KinematicInductionState<TGenTraits>::Velocity &velV)
   {
      SmartTruncation pTrunc = velV.oc().trunc();

      if(KinematicInductionState<TGenTraits>::VELOCITY_CASE == 20)
      {
         // Set toroidal component to zero
         for(int l=0; l < pTrunc->local()->spec()->nL(); ++l)
         {
            // Set some perturbation random energy in Toroidal component
            velV.rOc().rPerturbation().rTor().rLShell(l).setConstant(EPMComplex(0.0,0.0));
         }

         // **************************************************************
         // ******* WARNING THIS WILL NOT WORK IF DONE IN PARALLEL *******
         // **************************************************************
         // Set S20
         velV.rOc().rPerturbation().rPol().rLShell(2)(0,0) = EPMComplex(0.255832, 0.0);
         velV.rOc().rPerturbation().rPol().rLShell(2)(1,0) = EPMComplex(-0.2860290, 0.0);
      }

   }

   template <typename TGenTraits> KinematicInductionState<TGenTraits>::KinematicInductionState()
   {
   }
}

#endif // KINEMATICINDUCTIONSTATE_HPP
