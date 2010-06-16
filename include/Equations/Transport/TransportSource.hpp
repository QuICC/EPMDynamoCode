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
#include "Simulations/SimulationTraits.hpp"
#include "General/EPMTypedefs.hpp"
#include "Equations/Parameters/EquationParameters.hpp"
#include "Equations/Transport/TransportBase.hpp"

namespace EPMDynamo {

   /**
    * @brief General representation of the Transport equation
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class TransportSource: public TransportBase<TSimType, TSimTraits>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSim>::TransformType    TransformType;

         /**
         * @brief Constructor
         *
         * @param rC Codensity scalar (stored as reference)
         * @param transform Transform object (stored as reference)
         * \param tsteps Timestep parameters
         * @param params Simulation equation parameters
         */
         TransportSource(typename TSimTraits<TSim>::CodType &rC, TransformType &transform, TimestepParameters &tsteps, const EquationParameters &params);

         /**
         * @brief Simple empty destructor
         */
         virtual ~TransportSource() {};

         /**
          * @brief Update RTP of the equation
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

         /**
          * @brief Add codensity source term to the non linear part
          */
         void addSourceTerm();
         
      protected:

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits> TransportSource<TSimType, TSimTraits>::TransportSource(typename TSimTraits<TSim>::CodType &rC, typename TransportSource<TSimType, TSimTraits>::TransformType &transform, TimestepParameters &tsteps, const EquationParameters &params)
      : TransportBase<TSimType, TSimTraits>(rC, transform, tsteps, params)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void TransportSource<TSimType, TSimTraits>::addSourceTerm()
   {
      int nL = this->mrX.oc().perturbation().nL();
      const int l0 = this->mrX.oc().perturbation().minL();

      for(int l = l0; l < nL; ++l)
      {
         this->mNTerms.rOc().rPerturbation().rLShell(l) += this->mrX.src().lshell(l);
      }
   }

   template <typename TSimType, template <typename> class TSimTraits> void TransportSource<TSimType, TSimTraits>::updateRTP(const int step)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void TransportSource<TSimType, TSimTraits>::updateRHS()
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void TransportSource<TSimType, TSimTraits>::transformRHS(const int step)
   {
      // Transform Non Linear terms to spectral space from mNTerms values
      this->transformNTerms();
   }

}

#endif // TRANSPORTSOURCE_HPP
