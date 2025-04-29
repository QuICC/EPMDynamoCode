/** \file ThetaInfluenceCorrector.hpp
 *  \brief Implemenation of the corrector step of the theta method
 */

#ifndef THETAINFLUENCECORRECTOR_HPP
#define THETAINFLUENCECORRECTOR_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "Timestepping/PredictorCorrector/Theta/ThetaCorrector.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaLHSTOperatorSet.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaRHSTOperatorSet.hpp"
#include "Equations/InfluenceMatrix.hpp"

namespace EPMDynamo {

   /**
    * \brief Implemenation of the corrector step of the theta method
    */
   class ThetaInfluenceCorrector: public ThetaCorrector
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

         /// Typedef for a smart pointer to an influence matrix operator
         typedef EPMSHARED_PTR<InfluenceMatrix > SmartInfluenceMatrix;

         /**
          * @brief Constructor
          *
          * @param pOpLHS Pointer to the LHS operator
          * @param pOpRHS Pointer to the RHS operator
          */
         ThetaInfluenceCorrector(SmartInfluenceMatrix pInfluence, SmartScalarType pPrevious, SmartLHSOperators pOpLHS, SmartRHSOperators pOpRHS);

         /**
          * @brief Destructor
          */
         virtual ~ThetaInfluenceCorrector() {};

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
          * @brief Set the corrector step RHS
          *
          * @param newNTerms New non linear terms
          */
         void setRHS(ScalarType& rNTerms);

         /**
          * @brief Solve the corrector equations
          */
         void solve(ScalarType& rNTerms);

      private:
   };

}

#endif // THETAINFLUENCECORRECTOR_HPP
