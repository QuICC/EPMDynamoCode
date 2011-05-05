/** \file InductionLinear.hpp
 *  \brief Implementatio of the linear induction (with imposed field) equation
 */

#ifndef INDUCTIONLINEAR_HPP
#define INDUCTIONLINEAR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Config/SimulationConfig.hpp"
#include "General/EPMTypedefs.hpp"
#include "Equations/Induction/InductionBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementatio of the linear induction (with imposed field) equation
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the simulation implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class InductionLinear : public InductionBase<TSimType, TSimTraits>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef SimulationConfig::TransformType    TransformType;

         /// Typedef for the EquationParameters type
         typedef SimulationConfig::EquationParametersType EquationParametersType;

         /**
          * @brief Constructor
          *
          * @param rB Magnetic field (stored as reference)
          * @param rV Velocity field (stored as reference)
          * @param transform Transform object (stored as reference)
          * @param tsteps Timestep parameters
          * @param params Equation parameters
          */
         InductionLinear(typename TSimTraits<TSimType>::MagType &rB, typename TSimTraits<TSimType>::VelType &rV, TransformType &transform, TimestepParameters &tsteps, EquationParametersType &params);

         /**
          * @brief Simple empty destructor
          */
         virtual ~InductionLinear() {};

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

   template <typename TSimType, template <typename> class TSimTraits> InductionLinear<TSimType, TSimTraits>::InductionLinear(typename TSimTraits<TSimType>::MagType &rB, typename TSimTraits<TSimType>::VelType &rV, typename InductionLinear<TSimType, TSimTraits>::TransformType &transform, TimestepParameters &tsteps,  typename InductionLinear<TSimType, TSimTraits>::EquationParametersType &params)
      : InductionBase<TSimType, TSimTraits>(rB, transform, tsteps, params), mrV(rV)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void InductionLinear<TSimType, TSimTraits>::updateRTP(const int step)
   {
      // Update real space values for Magnetic field
      this->mrX.rOc().transform(step);

      // Update real space values for velocity field
      this->mrV.rOc().transform(step);
   }

   template <typename TSimType, template <typename> class TSimTraits> void InductionLinear<TSimType, TSimTraits>::updateRHS()
   {
      // Compute cross product \f$\vec{u}\times\vec{B}\f$
      this->mrV.oc().rtp().template cross<0>(this->mNTerms.rOc().rRTP(), this->mrX.oc().rtp(), this->mrParams.indAdvection());
   }

   template <typename TSimType, template <typename> class TSimTraits> void InductionLinear<TSimType, TSimTraits>::transformRHS(const int step)
   {
      // Transform Non Linear terms to spectral space from mNTerms values
      this->transformNTerms(this->mNTerms.rOc().rPerturbation().rPol(), this->mNTerms.rOc().rPerturbation().rTor());
   }
}

#endif // INDUCTIONLINEAR_HPP
