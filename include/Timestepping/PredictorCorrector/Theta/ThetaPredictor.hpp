/** \file ThetaPredictor.hpp
 *  \brief Implemenation of the predictor step of the theta method
 */

#ifndef THETAPREDICTOR_HPP
#define THETAPREDICTOR_HPP

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
#include "Timestepping/ErrorControl.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaLHSTOperatorSet.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaRHSTOperatorSet.hpp"

namespace EPMDynamo {

   /**
    * \brief Implemenation of the predictor step of the theta method
    */
   class ThetaPredictor: public SchemeIteration
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef SimulationConfig::NumericalScheme::ScalarType    ScalarType;

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
         ThetaPredictor(SmartLHSOperators pOpLHS, SmartRHSOperators pOpRHS);

         /**
          * @brief Destructor
          */
         virtual ~ThetaPredictor() {};

         /**
          * @brief Compute the implemented iteration
          *
          * @param rVar Input/Output variable
          * @param rNTerms New non linear terms
          */
         virtual void compute(ScalarType &rVar, ScalarType &rNTerms);
         
      protected:
         /**
          * @brief Pointer to the LHS operator
          */
         SmartLHSOperators  mpOpLHS;

         /**
          * @brief Pointer to the RHS operator
          */
         SmartRHSOperators  mpOpRHS;

         /**
          * @brief Set the predictor step RHS
          *
          * @param rVar Input/Output variable
          * @param oldVar Previous variable
          * @param nTerms Non linear terms
          */
         void setRHS(ScalarType& rVar, const ScalarType& oldVar, const ScalarType& nTerms);

         /**
          * @brief Solve the predictor equations
          *
          * @param rVar Input/Output variable
          */
         void solve(ScalarType& rVar);

      private:
   };

}

#endif // THETAPREDICTOR_HPP
