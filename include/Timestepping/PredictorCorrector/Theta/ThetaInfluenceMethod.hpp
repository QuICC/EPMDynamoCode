/** \file ThetaInfluenceMethod.hpp
 *  \brief Implementation of the predictor-corrector theta method with influence matrix
 */

#ifndef THETAINFLUENCEMETHOD_HPP
#define THETAINFLUENCEMETHOD_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "Timestepping/TimestepParameters.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaMethod.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaInfluencePredictor.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaInfluenceCorrector.hpp"
#include "Equations/InfluenceMatrix.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the predictor-corrector theta method with influence matrix
    */
   class ThetaInfluenceMethod: public ThetaMethod
   {
      public:
         /// Typedef from Simulation trait to local radial basis type
         typedef SimulationConfig::NumericalScheme::RadialBasisType    BasisType;

         /// Typedef from Simulation trait to local scalar type
         typedef SimulationConfig::NumericalScheme::ScalarType    ScalarType;

         /// Typedef for a smart pointer to an influence matrix operator
         typedef EPMSHARED_PTR<InfluenceMatrix > SmartInfluenceMatrix;

         /**
          * @brief Constructor
          *
          * @param a Coefficient \f$a\f$ of timestep scheme
          * @param b Coefficient \f$b\f$ of timestep scheme
          * @param basis Reference to the basis used for the operators
          * @param tsteps Timestep parameters
          * @param pTrunc Truncation information
          * @param hasL0 Is l=0 mode required?
          */
         ThetaInfluenceMethod(EPMFloat a, EPMFloat b, const BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0);

         /**
          * @brief Destructor
          */
         virtual ~ThetaInfluenceMethod() {};

         /**
          * @brief Add a boundary condition
          *
          * @param pBC Boundary condition
          */
         virtual void addBC(SmartBC pBC);

         /**
          * @brief Initialise the method
          */
         void init();
         
      protected:
         /**
          * @brief Update the timestep matrices after a timestep change
          */
         void updateTimeMatrices();

      private:
         /**
          * @brief Number of boundary conditions required for influence matrix
          */
         int mInfluenceNBC;

         /**
          * @brief Influence matrix
          */
         SmartInfluenceMatrix mpInfluence;

         /**
          * @brief Initialise the storage for the method
          */
         void initStorage();

         /**
          * @brief Initialise the influenc matrix
          */
         void initInfluence();

         /**
          * @brief update the influence matrix solution
          */
         void updateInfluence();
   };

}

#endif // THETAINFLUENCEMETHOD_HPP
