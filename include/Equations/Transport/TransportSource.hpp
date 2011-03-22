/** \file TransportSource.hpp
 *  \brief Implementation of transport diffusion with source term
 */

#ifndef TRANSPORTSOURCE_HPP
#define TRANSPORTSOURCE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "General/EPMTypedefs.hpp"
#include "Equations/Transport/TransportDiffusion.hpp"

namespace EPMDynamo {

   /**
    * @brief General representation of the Transport equation
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class TransportSource: public TransportDiffusion<TSimType, TSimTraits>
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
         TransportSource(typename TSimTraits<TSimType>::CodType &rC, TransformType &transform, TimestepParameters &tsteps, const EquationParametersType &params);

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

   template <typename TSimType, template <typename> class TSimTraits> TransportSource<TSimType, TSimTraits>::TransportSource(typename TSimTraits<TSimType>::CodType &rC, typename TransportSource<TSimType, TSimTraits>::TransformType &transform, TimestepParameters &tsteps, const  typename TransportSource<TSimType, TSimTraits>::EquationParametersType &params)
      : TransportDiffusion<TSimType, TSimTraits>(rC, transform, tsteps, params)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void TransportSource<TSimType, TSimTraits>::addSourceTerm()
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
