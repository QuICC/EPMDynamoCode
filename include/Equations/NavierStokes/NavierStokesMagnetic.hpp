/** \file NavierStokesMagnetic.hpp
 *  \brief Implementation of the Magneto-convection Navier-Stokes equation
 */

#ifndef NAVIERSTOKESMAGNETIC_HPP
#define NAVIERSTOKESMAGNETIC_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/SimulationTraits.hpp"
#include "General/EPMTypedefs.hpp"
#include "Equations/NavierStokes/NavierStokesRotating.hpp"
#include "Equations/EquationParameters.hpp"

namespace EPMDynamo {

   /**
    * @brief This class implements the Navier-Stokes equation including the magnetic effects
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the simulation implementation
    * \tparam TTStepPol NEED TO CHANGE THIS
    *
    * \bug Need to change the implementation of the poisson solver
    */
   template <typename TSimType, template <typename> class TSimTraits> class NavierStokesMagnetic : public NavierStokesRotating<TSimType, TSimTraits>  {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSim>::TransformType    TransformType;

         /// Typedef from Simulation trait to local truncation type
         typedef typename TSim::ScalarType    ScalarType;

         /**
          * @brief Constructor
          *
          * @param rV Magnetic field (stored as reference)
          * @param rB Velocity field (stored as reference)
          * @param transform Transform object (stored as reference)
          * \param tsteps Timestep parameters
          * @param params Simulation equation paramters
          */
         NavierStokesMagnetic(typename TSimTraits<TSim>::VelType &rV, typename TSimTraits<TSim>::MagType &rB, TransformType &transform, TimestepParameters &tsteps, EquationParameters &params);

         /**
          * @brief Simple empty destructor
          */
         virtual ~NavierStokesMagnetic();

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
         typename TSimTraits<TSim>::MagType&  mrB;

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits> NavierStokesMagnetic<TSimType, TSimTraits>::NavierStokesMagnetic(typename TSimTraits<TSim>::VelType &rV, typename TSimTraits<TSim>::MagType &rB, typename NavierStokesMagnetic<TSimType, TSimTraits>::TransformType &transform, TimestepParameters &tsteps, EquationParameters &params)
      : NavierStokesRotating<TSimType, TSimTraits>(rV, transform, tsteps, params), mrB(rB)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesMagnetic<TSimType, TSimTraits>::updateRTP(const int step)
   {
      // Update real space values of velocity field
      this->mrX.rOc().transform(step);

      // Update real space values of the curl of the velocity field
      this->mrX.rOc().curlTransform(step);

      // Update real space values of magnetic field
      this->mrB.rOc().transform(step);

      // Update real space values of the curl of the magnetic field
      this->mrB.rOc().curlTransform(step);
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesMagnetic<TSimType, TSimTraits>::updateRHS()
   {
      // Compute \f$u\times (\nabla \times u) \f$
      this->mrX.oc().rtp().template cross<0>(this->mNTerms.rOc().rRTP(), this->mrX.oc().curl(), this->mrParams.Ro());

      // Compute \f$(\nabla \times B)\times B\f$
      this->mrB.oc().curl().template cross<1>(this->mNTerms.rOc().rRTP(), this->mrB.oc().rtp());

      // Compute \f$\hat{z}\times\vec{u}\f$
      this->mrX.oc().rtp().template crossZVect<-1>(this->mNTerms.rOc().rRTP());
   }

}

#endif // NAVIERSTOKESMAGNETIC_HPP
