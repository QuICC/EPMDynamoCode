/** \file NavierStokesMHD.hpp
 *  \brief Implementation of the MHD Navier-Stokes equation
 */

#ifndef NAVIERSTOKESMHD_HPP
#define NAVIERSTOKESMHD_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "General/EPMTypedefs.hpp"
#include "Equations/NavierStokes/NavierStokesMagnetic.hpp"
#include "Equations/Parameters/EquationParameters.hpp"

namespace EPMDynamo {

   /**
    * @brief This class implements the full MHD Navierr-Stokes equation
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class NavierStokesMHD : public NavierStokesMagnetic<TSimType, TSimTraits>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::ScalarType    ScalarType;

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
         NavierStokesMHD(typename TSimTraits<TSimType>::VelType &rV, typename TSimTraits<TSimType>::MagType &rB, typename TSimTraits<TSimType>::CodType &rC, TransformType &transform, TimestepParameters &tsteps, EquationParameters &params);

         /**
          * @brief Simple empty destructor
          */
         virtual ~NavierStokesMHD() {};

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
          * @brief Const Reference variable to the codensity scalar
          */
         typename TSimTraits<TSimType>::CodType&  mrC;

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits> NavierStokesMHD<TSimType, TSimTraits>::NavierStokesMHD(typename TSimTraits<TSimType>::VelType &rV, typename TSimTraits<TSimType>::MagType &rB, typename TSimTraits<TSimType>::CodType &rC, typename NavierStokesMHD<TSimType, TSimTraits>::TransformType &transform, TimestepParameters &tsteps, EquationParameters &params)
      : NavierStokesMagnetic<TSimType, TSimTraits>(rV, rB, transform, tsteps, params), mrC(rC)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesMHD<TSimType, TSimTraits>::updateRTP(const int step)
   {
      // Update real space values of velocity field
      this->mrX.rOc().transform(step);

      // Update real space values of the curl of the velocity field
      this->mrX.rOc().curlTransform(step);

      // Update real space values of magnetic field
      this->mrB.rOc().transform(step);

      // Update real space values of the curl of the magnetic field
      this->mrB.rOc().curlTransform(step);

      // Update real space values of codensity scalar
      this->mrC.rOc().transform(step);
   }

   template <typename TSimType, template <typename> class TSimTraits> void NavierStokesMHD<TSimType, TSimTraits>::updateRHS()
   {
      // Compute \f$u\times (\nabla \times B) \f$
      this->mrX.oc().rtp().template cross<0>(this->mNTerms.rOc().rRTP(), this->mrX.oc().curl(), this->mrParams.Ro());

      // Compute \f$(\nabla \times B)\times B\f$
      this->mrB.oc().curl().template cross<1>(this->mNTerms.rOc().rRTP(), this->mrB.oc().rtp());

      // Compute \f$C \vec{r}\f$
      this->mrC.oc().rtp().template radVect<1>(this->mNTerms.rOc().rRTP(), this->mrParams.q()*this->mrParams.Ra());

      // Compute \f$\hat{z}\times\vec{u}\f$
      this->mrX.oc().rtp().template crossZVect<-1>(this->mNTerms.rOc().rRTP());
   }

}

#endif // NAVIERSTOKESMHD_HPP
