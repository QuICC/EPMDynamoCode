/** \file TransportDummy.hpp
 *  \brief Implementation of the DUMMY transport equation
 */

#ifndef TRANSPORTDUMMY_HPP
#define TRANSPORTDUMMY_HPP

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
    * @brief This class implements the full DUMMY Transport equation
    *
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimTraits> class TransportDummy: public TransportSource<TSimTraits>
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
         TransportDummy(typename TSimTraits::CodType &rC, typename TSimTraits::VelType &rV, TransformType &transform, TimestepParameters &tsteps, const EquationParametersType &params);

         /**
         * @brief Simple empty destructor
         */
         virtual ~TransportDummy() {};

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

   template <typename TSimTraits> TransportDummy<TSimTraits>::TransportDummy(typename TSimTraits::CodType &rC, typename TSimTraits::VelType &rV, typename TransportDummy<TSimTraits>::TransformType &transform, TimestepParameters &tsteps, const  typename TransportDummy<TSimTraits>::EquationParametersType &params)
      : TransportSource<TSimTraits>(rC, transform, tsteps, params), mrV(rV)
   {
   }

   template <typename TSimTraits> void TransportDummy<TSimTraits>::updateRTP(const int step)
   {
      // Update real space values of the gradient of codensity scalar
      this->mrX.rOc().gradTransform(step);

      // Update real space values for velocity field
      this->mrV.rOc().transform(step);
   }

   template <typename TSimTraits> void TransportDummy<TSimTraits>::updateRHS()
   {
      // Compute \f$u\cdot\nabla C\f$
      this->mrV.oc().rtp().template dot<0>(this->mNTerms.rOc().rRTP(), this->mrX.oc().grad(), -this->mrParams.tptAdvection());
   }

   template <typename TSimTraits> void TransportDummy<TSimTraits>::timestep()
   {
      this->mrX.rOc().rPerturbation();
   }

}

#endif // TRANSPORTDUMMY_HPP
