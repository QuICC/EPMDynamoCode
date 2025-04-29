/** \file ThetaCorrector.hpp
 *  \brief Implemenation of the corrector step of the theta method
 */

#ifndef THETACORRECTOR_HPP
#define THETACORRECTOR_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "Timestepping/SchemeIteration.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaLHSTOperatorSet.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaRHSTOperatorSet.hpp"

namespace EPMDynamo {

   /**
    * \brief Implemenation of the corrector step of the theta method
    */
   class ThetaCorrector: public SchemeIteration
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef SimulationConfig::NumericalScheme::ScalarType    ScalarType;

         /// typedef for a pointer to a scalar type
          typedef EPMSHARED_PTR<ScalarType> SmartScalarType;

         /// Typedef for a smart pointer to LHS operator set
         typedef EPMSHARED_PTR<ThetaLHSTOperatorSet > SmartLHSOperators;

         /// Typedef for a smart pointer to RHS operator set
         typedef EPMSHARED_PTR<ThetaRHSTOperatorSet > SmartRHSOperators;

         /**
          * @brief Constructor
          *
          * @param pOpLHS Pointer to the LHS operator
          * @param pOpRHS Pointer to the RHS operator
          */
         ThetaCorrector(SmartScalarType pPrevious, SmartLHSOperators pOpLHS, SmartRHSOperators pOpRHS);

         /**
          * @brief Destructor
          */
         virtual ~ThetaCorrector() {};

         /**
          * @brief Compute the implemented iteration
          *
          * @param rVar Input/Output variable
          * @param rNTerms New non linear terms
          */
         virtual void compute(ScalarType &rVar, ScalarType &rNTerms);

         /**
          * @brief Does this iteration provide timestep error information
          */
         bool providesError() const;

         /**
          * @brief Keep the system in the corrector iteration
          */
         bool doNextIteration() const;
         
      protected:
         /**
          * @brief Pointer to the previous NTerms scalar value
          */
         SmartScalarType   mpPreviousNTerms;

         /**
          * @brief Pointer to the LHS operator
          */
         SmartLHSOperators  mpOpLHS;

         /**
          * @brief Pointer to the RHS operator
          */
         SmartRHSOperators  mpOpRHS;

         /**
          * @brief Set the corrector step RHS
          *
          * @param newNTerms New non linear terms
          */
         void setRHS(ScalarType& rNTerms);

         /**
          * @brief Solve the corrector equations
          */
         void solve(ScalarType& rNTerms);

         /**
          * @brief Use the corrector solution
          *
          * @param rVar Input/Output variable
          */
         void useCorrection(ScalarType& rVar, const ScalarType& nTerms);

         /**
          * @brief Add correction from corrector step
          *
          * @param rVar Input/Output variable
          * @param corr Corrector correction
          */
         void addCorrection(ScalarType& rVar, const ScalarType& corr) const;

      private:
   };

   inline bool ThetaCorrector::providesError() const
   {
      return true;
   }

   inline bool ThetaCorrector::doNextIteration() const
   {
      return false;
   }

}

#endif // THETACORRECTOR_HPP
