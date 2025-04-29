/** \file SimulationBase.hpp
 *  \brief Implementation independent parts of a simulation
 */

#ifndef SIMULATIONBASE_HPP
#define SIMULATIONBASE_HPP

// Configuration includes
//
#include "Config/Parallelisation.h"
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Simulations/ComputationBase.hpp"
#include "Simulations/SimulationControl.hpp"
#include "IO/ASCII/FDSHSplittingFile.hpp"
#include "IO/ASCII/PhysicalSplittingFile.hpp"
#include "IO/ASCII/SpectralSplittingFile.hpp"
#include "IO/IOSystem.hpp"

namespace EPMDynamo {

   /**
    * @brief Implements the implementation independent part of a simulation
    */
   class SimulationBase: public ComputationBase<IOSystem>
   {
      public:
         /// Typedef for the transform type
         typedef SimulationConfig::TransformType  TransformType;

         /// Typedef for the EquationParameters type
         typedef SimulationConfig::EquationParametersType EquationParametersType;

         /**
          * @brief Simple empty destructor
          */
         virtual ~SimulationBase() {};

      protected:
         /**
          * @brief Constructor
          */
         SimulationBase();

         /**
          * @brief SimulationBase control
          */
         SimulationControl  mSimControl;

         /**
          * @brief write data to files
          *
          * \epmTodo Rename this to a more explicit name and review the used approach
          */
         void writeFiles();

         /**
          * @brief Initialise the IO output
          */
         void initOutput();

         /**
          * @brief Do operations required just before starting the time integration
          */
         void preRun();

         /**
          * @brief Do operations required just after finishing the time integration
          */
         void postRun();

         /**
         * @name Methods that need to be overloaded in implementation
         */
         //@{

         /**
          * @brief Initialise the fields
          */
         void initFields();

         /**
          * @brief Add HDF5 output files
          */
         void addHDF5Output();

         /**
          * @brief Add ASCII output files
          */
         void addASCIIOutput();

         /**
          * @brief Initialise the equations
          */
         void initEquations();

         /**
          * @brief Update equations RTP values
          *
          * \param step Step in the transform for multistep trransforms
          */
         void updateEquationsRTP(const int step);

         /**
          * @brief Update equations RHS
          */
         void updateEquationsRHS();

         /**
          * @brief Transform equations RHS
          *
          * \param step Step in the transform for multistep trransforms
          */
         void transformEquationsRHS(const int step);

         /**
          * @brief Include external influence
          */
         void addExternalInfluence();

         /**
          * @brief Timestep the equations
          */
         void timestepEquations();
         //@}
         
         /**
          * @brief Create Splitting visualisation files
          */
         void describeSplitting();

      private:
   };
}

#endif // SIMULATIONBASE_HPP
