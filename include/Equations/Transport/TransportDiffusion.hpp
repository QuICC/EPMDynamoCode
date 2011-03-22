/** \file TransportDiffusion.hpp
 *  \brief Implementation of the transport diffusion equation
 */

#ifndef TRANSPORTDIFFUSION_HPP
#define TRANSPORTDIFFUSION_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "General/EPMTypedefs.hpp"
#include "Equations/Transport/TransportBase.hpp"

namespace EPMDynamo {

   /**
    * @brief General representation of the Transport equation with only diffusion included
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class TransportDiffusion: public TransportBase<TSimType, TSimTraits>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /// Typedef for the EquationParameters type
         typedef typename SimulationTraits<TSimType>::EquationParametersType EquationParametersType;

         /**
         * @brief Constructor
         *
         * @param rC Codensity scalar (stored as reference)
         * @param transform Transform object (stored as reference)
         * \param tsteps Timestep parameters
         * @param params Simulation equation parameters
         */
         TransportDiffusion(typename TSimTraits<TSimType>::CodType &rC, TransformType &transform, TimestepParameters &tsteps, const EquationParametersType &params);

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

   template <typename TSimType, template <typename> class TSimTraits> TransportDiffusion<TSimType, TSimTraits>::TransportDiffusion(typename TSimTraits<TSimType>::CodType &rC, typename TransportDiffusion<TSimType, TSimTraits>::TransformType &transform, TimestepParameters &tsteps, const  typename TransportDiffusion<TSimType, TSimTraits>::EquationParametersType &params)
      : TransportBase<TSimType, TSimTraits>(rC, transform, tsteps, params)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void TransportDiffusion<TSimType, TSimTraits>::updateRTP(const int step)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void TransportDiffusion<TSimType, TSimTraits>::updateRHS()
   {
      // Set the RTP non linear terms to zero
      this->mNTerms.rOc().rRTP().initialiseZeros();
   }

   template <typename TSimType, template <typename> class TSimTraits> void TransportDiffusion<TSimType, TSimTraits>::transformRHS(const int step)
   {
      // Transform Non Linear terms to spectral space from mNTerms values
      this->transformNTerms();
   }

}

#endif // TRANSPORTDIFFUSION_HPP
