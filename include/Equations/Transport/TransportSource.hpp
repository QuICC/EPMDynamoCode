/** \file TransportSource.hpp
 *  \brief Implementation of transport diffusion with source term
 */

#ifndef TRANSPORTSOURCE_HPP
#define TRANSPORTSOURCE_HPP

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
#include "Equations/Transport/TransportDiffusion.hpp"

namespace EPMDynamo {

   /**
    * @brief General representation of the Transport equation
    *
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimTraits> class TransportSource: public TransportDiffusion<TSimTraits>
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
         TransportSource(typename TSimTraits::CodType &rC, TransformType &transform, TimestepParameters &tsteps, const EquationParametersType &params);

         /**
         * @brief Simple empty destructor
         */
         virtual ~TransportSource() {};

         /**
          * @brief Add codensity source term to the non linear part
          */
         void addSourceTerm();
         
      protected:

      private:
   };

   template <typename TSimTraits> TransportSource<TSimTraits>::TransportSource(typename TSimTraits::CodType &rC, typename TransportSource<TSimTraits>::TransformType &transform, TimestepParameters &tsteps, const  typename TransportSource<TSimTraits>::EquationParametersType &params)
      : TransportDiffusion<TSimTraits>(rC, transform, tsteps, params)
   {
   }

   template <typename TSimTraits> void TransportSource<TSimTraits>::addSourceTerm()
   {
      int nL = this->mrX.oc().perturbation().nL();
      const int l0 = this->mrX.oc().perturbation().minL();

      for(int l = l0; l < nL; ++l)
      {
         this->mNTerms.rOc().rPerturbation().rLShell(l) += this->mrX.ocSrc().lshell(l);
      }
   }
}

#endif // TRANSPORTSOURCE_HPP
