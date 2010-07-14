/** \file TransportLinear.hpp
 *  \brief Implementation of the linear transport equation
 */

#ifndef TRANSPORTLINEAR_HPP
#define TRANSPORTLINEAR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "General/EPMTypedefs.hpp"
#include "Equations/Parameters/EquationParameters.hpp"
#include "Equations/Transport/TransportSource.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the linear transport equation
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class TransportLinear: public TransportSource<TSimType, TSimTraits>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /**
         * @brief Constructor
         *
         * @param rC Codensity scalar (stored as reference)
         * @param rV Velocity field (stored as reference)
         * @param transform Transform object (stored as reference)
         * \param tsteps Timestep parameters
         * @param params Simulation equation parameters
         */
         TransportLinear(typename TSimTraits<TSimType>::CodType &rC, typename TSimTraits<TSimType>::VelType &rV, TransformType &transform, TimestepParameters &tsteps, const EquationParameters &params);

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
         typename TSimTraits<TSimType>::VelType&  mrV;

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits> TransportLinear<TSimType, TSimTraits>::TransportLinear(typename TSimTraits<TSimType>::CodType &rC, typename TSimTraits<TSimType>::VelType &rV, typename TransportLinear<TSimType, TSimTraits>::TransformType &transform, TimestepParameters &tsteps, const EquationParameters &params)
      : TransportSource<TSimType, TSimTraits>(rC, transform, tsteps, params), mrV(rV)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void TransportLinear<TSimType, TSimTraits>::updateRTP(const int step)
   {
      // Update real space values for velocity field
      this->mrV.rOc().transform(step);
   }

   template <typename TSimType, template <typename> class TSimTraits> void TransportLinear<TSimType, TSimTraits>::updateRHS()
   {
      // Compute \f$u\cdot\hat{r}\f$
      this->mrV.oc().rtp().r().template radVectProj<0>(this->mNTerms.rOc().rRTP(), -1.0);
   }

}

#endif // TRANSPORTLINEAR_HPP
