/** \file InductionMHD.hpp
 *  \brief Implementatio of the MHD induction equation
 */

#ifndef INDUCTIONMHD_HPP
#define INDUCTIONMHD_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "General/EPMTypedefs.hpp"
#include "Equations/Induction/InductionBase.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of the full MHD Induction equation
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class InductionMHD : public InductionBase<TSimType, TSimTraits>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /**
          * @brief Constructor
          *
          * @param rB Magnetic field (stored as reference)
          * @param rV Velocity field (stored as reference)
          * @param transform Transform object (stored as reference)
          * \param tsteps Timestep parameters
          */
         InductionMHD(typename TSimTraits<TSimType>::MagType &rB, typename TSimTraits<TSimType>::VelType &rV, TransformType &transform, TimestepParameters &tsteps);

         /**
          * @brief Simple empty destructor
          */
         virtual ~InductionMHD() {};

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

         /**
          * @brief Const Reference variable to the velocity field
          */
         typename TSimTraits<TSimType>::VelType&  mrV;

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits> InductionMHD<TSimType, TSimTraits>::InductionMHD(typename TSimTraits<TSimType>::MagType &rB, typename TSimTraits<TSimType>::VelType &rV, typename InductionMHD<TSimType, TSimTraits>::TransformType &transform, TimestepParameters &tsteps)
      : InductionBase<TSimType, TSimTraits>(rB, transform, tsteps), mrV(rV)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void InductionMHD<TSimType, TSimTraits>::updateRTP(const int step)
   {
      // Update real space values for Magnetic field
      this->mrX.rOc().transform(step);

      // Update real space values for velocity field
      this->mrV.rOc().transform(step);
   }

   template <typename TSimType, template <typename> class TSimTraits> void InductionMHD<TSimType, TSimTraits>::updateRHS()
   {
      // Compute cross product \f$\vec{u}\times\vec{B}\f$
      this->mrV.oc().rtp().template cross<0>(this->mNTerms.rOc().rRTP(), this->mrX.oc().rtp());
   }

   template <typename TSimType, template <typename> class TSimTraits> void InductionMHD<TSimType, TSimTraits>::transformRHS(const int step)
   {
      // Transform Non Linear terms to spectral space from mNTerms values
      this->transformNTerms(this->mNTerms.rOc().rPerturbation().rPol(), this->mNTerms.rOc().rPerturbation().rTor());
   }
}

#endif // INDUCTIONMHD_HPP
