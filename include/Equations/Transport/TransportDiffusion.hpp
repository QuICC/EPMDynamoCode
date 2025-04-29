/** \file TransportDiffusion.hpp
 *  \brief Implementation of the transport diffusion equation
 */

#ifndef TRANSPORTDIFFUSION_HPP
#define TRANSPORTDIFFUSION_HPP

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
#include "Equations/Transport/TransportBase.hpp"

namespace EPMDynamo {

   /**
    * @brief General representation of the Transport equation with only diffusion included
    *
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimTraits> class TransportDiffusion: public TransportBase<TSimTraits>
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
         * @param transform Transform object (stored as reference)
         * \param tsteps Timestep parameters
         * @param params Simulation equation parameters
         */
         TransportDiffusion(typename TSimTraits::CodType &rC, TransformType &transform, TimestepParameters &tsteps, const EquationParametersType &params);

         /**
         * @brief Simple empty destructor
         */
         virtual ~TransportDiffusion() {};

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

   template <typename TSimTraits> TransportDiffusion<TSimTraits>::TransportDiffusion(typename TSimTraits::CodType &rC, typename TransportDiffusion<TSimTraits>::TransformType &transform, TimestepParameters &tsteps, const  typename TransportDiffusion<TSimTraits>::EquationParametersType &params)
      : TransportBase<TSimTraits>(rC, transform, tsteps, params)
   {
   }

   template <typename TSimTraits> void TransportDiffusion<TSimTraits>::updateRTP(const int step)
   {
   }

   template <typename TSimTraits> void TransportDiffusion<TSimTraits>::updateRHS()
   {
      // Set the RTP non linear terms to zero
      this->mNTerms.rOc().rRTP().initialiseZeros();
   }

   template <typename TSimTraits> void TransportDiffusion<TSimTraits>::transformRHS(const int step)
   {
      // Transform Non Linear terms to spectral space from mNTerms values
      this->transformNTerms();
   }

}

#endif // TRANSPORTDIFFUSION_HPP
