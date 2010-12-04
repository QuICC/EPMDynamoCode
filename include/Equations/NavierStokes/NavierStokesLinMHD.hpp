/** \file NavierStokesLinMHD.hpp
 *  \brief Implementation of the linear MHD Navier-Stokes equation
 */

#ifndef NAVIERSTOKESLINMHD_HPP
#define NAVIERSTOKESLINMHD_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "General/EPMTypedefs.hpp"
#include "Equations/NavierStokes/NavierStokesRotating.hpp"

namespace EPMDynamo {

   /**
    * @brief This class implements the linear MHD Navierr-Stokes equation
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class NavierStokesLinMHD : public NavierStokesRotating<TSimType, TSimTraits>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::ScalarType    ScalarType;

         /// Typedef for the EquationParameters type
         typedef typename SimulationTraits<TSimType>::EquationParametersType EquationParametersType;

         /**
          * @brief Constructor
          *
          * @param rV Magnetic field (stored as reference)
          * @param rB Velocity field (stored as reference)
          * @param rC Codensity scalar (stored as reference)
          * @param transform Transform object (stored as reference)
          * \param tsteps Timestep parameters
          * @param params Simulation equation paramters
          */
         NavierStokesLinMHD(typename TSimTraits<TSimType>::VelType &rV, typename TSimTraits<TSimType>::MagType &rB, typename TSimTraits<TSimType>::CodType &rC, TransformType &transform, TimestepParameters &tsteps, EquationParametersType &params);

         /**
          * @brief Simple empty destructor
          */
         virtual ~NavierStokesLinMHD() {};

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
         
      protected:

         /**
          * @brief Const Reference variable to the magnetic field
          */
         typename TSimTraits<TSimType>::MagType&  mrB;

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits> NavierStokesLinMHD<TSimType, TSimTraits>::NavierStokesLinMHD(typename TSimTraits<TSimType>::VelType &rV, typename TSimTraits<TSimType>::MagType &rB, typename TSimTraits<TSimType>::CodType &rC, typename NavierStokesLinMHD<TSimType, TSimTraits>::TransformType &transform, TimestepParameters &tsteps,  typename NavierStokesLinMHD<TSimType, TSimTraits>::EquationParametersType &params)
      : NavierStokesRotating<TSimType, TSimTraits>(rV, rC, transform, tsteps, params), mrB(rB)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesLinMHD<TSimType, TSimTraits>::updateRTP(const int step)
   {
      // Update real space values of velocity field
      this->mrX.rOc().transform(step);

      // Update real space values of magnetic field
      this->mrB.rOc().transform(step);

      // Update real space values of codensity scalar
      this->mrC.rOc().transform(step);
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesLinMHD<TSimType, TSimTraits>::updateRHS()
   {
      // Compute \f$(\nabla \times B)\times B\f$
      this->mrB.oc().curl().template cross<1>(this->mNTerms.rOc().rRTP(), this->mrB.oc().rtp(), this->mrParams.nsLorentz());

      // Compute \f$C \vec{r}\f$
      this->mrC.oc().rtp().template radVect<1>(this->mNTerms.rOc().rRTP(), this->mrParams.nsBuoyancy());

      // Compute \f$\hat{z}\times\vec{u}\f$
      this->mrX.oc().rtp().template crossZVect<-1>(this->mNTerms.rOc().rRTP(), this->mrParams.nsCoriolis());
   }

}

#endif // NAVIERSTOKESLINMHD_HPP
