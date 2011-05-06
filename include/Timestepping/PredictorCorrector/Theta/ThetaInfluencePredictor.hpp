/** \file ThetaInfluencePredictor.hpp
 *  \brief Implemenation of the predictor step of the theta method
 */

#ifndef THETAINFLUENCePREDICTOR_HPP
#define THETAINFLUENCEPREDICTOR_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "Timestepping/PredictorCorrector/Theta/ThetaPredictor.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaLHSTOperatorSet.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaRHSTOperatorSet.hpp"
#include "Equations/InfluenceMatrix.hpp"

namespace EPMDynamo {

   /**
    * \brief Implemenation of the predictor step of the theta method
    */
   class ThetaInfluencePredictor: public ThetaPredictor
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef SimulationConfig::NumericalScheme::ScalarType    ScalarType;

         /// Typedef for a smart pointer to LHS operator set
         typedef EPMSHARED_PTR<ThetaLHSTOperatorSet > SmartLHSOperators;

         /// Typedef for a smart pointer to RHS operator set
         typedef EPMSHARED_PTR<ThetaRHSTOperatorSet > SmartRHSOperators;

         /// Typedef for a smart pointer to an influence matrix operator
         typedef EPMSHARED_PTR<InfluenceMatrix > SmartInfluenceMatrix;

         /**
          * @brief Constructor
          *
          * @param pOpLHS Pointer to the LHS operator
          * @param pOpRHS Pointer to the RHS operator
          */
         ThetaInfluencePredictor(SmartInfluenceMatrix pInfluence, SmartLHSOperators pOpLHS, SmartRHSOperators pOpRHS);

         /**
          * @brief Destructor
          */
         virtual ~ThetaInfluencePredictor() {};

         /**
          * @brief Compute the implemented iteration
          *
          * @param rVar Input/Output variable
          * @param rNTerms New non linear terms
          */
         virtual void compute(ScalarType &rVar, ScalarType &rNTerms);
         
      protected:
         /**
          * @brief Pointer to the influence matrix operator
          */
         SmartInfluenceMatrix  mpInfluence;

         /**
          * @brief Set the predictor step RHS
          *
          * @param rVar Input/Output variable
          * @param oldVar Previous variable
          * @param nTerms Non linear terms
          */
         void setRHS(ScalarType& rVar, const ScalarType& oldVar, ScalarType& rNTerms);

         /**
          * @brief Solve the predictor equations
          *
          * @param rVar Input/Output variable
          */
         void solve(ScalarType& rVar);

      private:
   };

}

#endif // THETAINFLUENCEPREDICTOR_HPP
