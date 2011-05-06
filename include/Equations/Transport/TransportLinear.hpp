/** \file TransportLinear.hpp
 *  \brief Implementation of the linear transport equation
 */

#ifndef TRANSPORTLINEAR_HPP
#define TRANSPORTLINEAR_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Equations/Transport/TransportSource.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the linear transport equation
    *
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimTraits> class TransportLinear: public TransportSource<TSimTraits>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef SimulationConfig::TransformType    TransformType;

         /// Typedef for the EquationParameters type
         typedef SimulationConfig::EquationParametersType EquationParametersType;

         /**
         * @brief Constructor
         *
         * @param rC Codensity scalar (stored as reference)
         * @param rV Velocity field (stored as reference)
         * @param transform Transform object (stored as reference)
         * \param tsteps Timestep parameters
         * @param params Simulation equation parameters
         */
         TransportLinear(typename TSimTraits::CodType &rC, typename TSimTraits::VelType &rV, TransformType &transform, TimestepParameters &tsteps, const EquationParametersType &params);

         /**
         * @brief Simple empty destructor
         */
         virtual ~TransportLinear() {};

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
          * @brief Const Reference variable to the velocity field
          */
         typename TSimTraits::VelType&  mrV;

      private:
   };

   template <typename TSimTraits> TransportLinear<TSimTraits>::TransportLinear(typename TSimTraits::CodType &rC, typename TSimTraits::VelType &rV, typename TransportLinear<TSimTraits>::TransformType &transform, TimestepParameters &tsteps, const  typename TransportLinear<TSimTraits>::EquationParametersType &params)
      : TransportSource<TSimTraits>(rC, transform, tsteps, params), mrV(rV)
   {
   }

   template <typename TSimTraits> void TransportLinear<TSimTraits>::updateRTP(const int step)
   {
      // Update real space values for velocity field
      this->mrV.rOc().transform(step);
   }

   template <typename TSimTraits> void TransportLinear<TSimTraits>::updateRHS()
   {
      // Compute \f$u\cdot\hat{r}\f$
      this->mrV.oc().rtp().r().template radVectProj<0>(this->mNTerms.rOc().rRTP(), -this->mrParams.tptAdvection());
   }

}

#endif // TRANSPORTLINEAR_HPP
