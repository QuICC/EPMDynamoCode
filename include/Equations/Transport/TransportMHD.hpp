/** \file TransportMHD.hpp
 *  \brief Implementation of the MHD transport equation
 */

#ifndef TRANSPORTMHD_HPP
#define TRANSPORTMHD_HPP

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
    * @brief This class implements the full MHD Transport equation
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class TransportMHD: public TransportSource<TSimType, TSimTraits>
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
         TransportMHD(typename TSimTraits<TSimType>::CodType &rC, typename TSimTraits<TSimType>::VelType &rV, TransformType &transform, TimestepParameters &tsteps, const EquationParameters &params);

         /**
         * @brief Simple empty destructor
         */
         virtual ~TransportMHD() {};

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

   template <typename TSimType, template <typename> class TSimTraits> TransportMHD<TSimType, TSimTraits>::TransportMHD(typename TSimTraits<TSimType>::CodType &rC, typename TSimTraits<TSimType>::VelType &rV, typename TransportMHD<TSimType, TSimTraits>::TransformType &transform, TimestepParameters &tsteps, const EquationParameters &params)
      : TransportSource<TSimType, TSimTraits>(rC, transform, tsteps, params), mrV(rV)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void TransportMHD<TSimType, TSimTraits>::updateRTP(const int step)
   {
      // Update real space values of the gradient of codensity scalar
      this->mrX.rOc().gradTransform(step);

      // Update real space values for velocity field
      this->mrV.rOc().transform(step);
   }

   template <typename TSimType, template <typename> class TSimTraits> void TransportMHD<TSimType, TSimTraits>::updateRHS()
   {
      // Compute \f$u\cdot\nabla C\f$
      this->mrV.oc().rtp().template dot<0>(this->mNTerms.rOc().rRTP(), this->mrX.oc().grad(), -1.0);
   }

}

#endif // TRANSPORTMHD_HPP
