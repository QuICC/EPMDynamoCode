/** \file SimulationBase.hpp
 *  \brief Implementation independent parts of a simulation
 */

#ifndef SIMULATIONBASE_HPP
#define SIMULATIONBASE_HPP

// Configuration includes
//
#include "Config/Parallelisation.h"

// System includes
//
#include <boost/static_assert.hpp>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Simulations/ComputationBase.hpp"
#include "Simulations/SimulationControl.hpp"
#include "IO/IOSystem.hpp"

namespace EPMDynamo {

   /**
    * @brief Implements the implementation independent part of a simulation
    *
    * \tparam TSimType General simulation type
    */
   template <typename TSimType> class SimulationBase: public ComputationBase<TSimType, IOSystem>
   {
      public:
         /// Typedef for the transform type
         typedef  typename SimulationTraits<TSimType>::TransformType  TransformType;

         /// Typedef for the EquationParameters type
         typedef typename SimulationTraits<TSimType>::EquationParametersType EquationParametersType;

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
         SimulationControl<TSimType>   mSimControl;

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

      private:
   };

   template <typename TSimType> SimulationBase<TSimType>::SimulationBase()
      : ComputationBase<TSimType, IOSystem>(), mSimControl(this->mIOSys.cfg()->aTStep(), this->mEqParams, this->mIOSys.cfg()->aRunI(), this->mIOSys.cfg()->aRun())
   {
   }

   template <typename TSimType> void SimulationBase<TSimType>::initOutput()
   {
      // Initialise all the create writers
      this->mIOSys.initWriters();
   }

   template <typename TSimType> void SimulationBase<TSimType>::writeFiles()
   {
      if(this->mSimControl.tsParams().isNextStep())
      {
         // Save ASCII files
         if(this->mSimControl.tsCounter().triggerASCII())
         {
            this->mIOSys.writeASCII();

            this->mSimControl.checkInterface();
         }

         // Save state file
         if(this->mSimControl.tsCounter().triggerState())
         {
            this->mIOSys.writeHDF5();
         }
      }
   }

   template <typename TSimType> void SimulationBase<TSimType>::preRun()
   {
      // Write loaded state to HDF5 file (this can be important if the truncations were not the same)
      this->mIOSys.writeHDF5();

      // Synchronize simulation
      EPMDYNAMO_SYNCHRONIZE;

      // Execute last minute simulation control initialisation operations
      this->mSimControl.preRun();
   }

   template <typename TSimType> void SimulationBase<TSimType>::postRun()
   {
      // Write the state to HDF5 file before finishing
      this->mIOSys.writeHDF5();

      // Execute early simulation control finalisation operations
      this->mSimControl.postRun();
   }

   template <typename TSimType> void SimulationBase<TSimType>::initFields()
   {
      BOOST_STATIC_ASSERT(sizeof(TSimType) == 0); 
   }

   template <typename TSimType> void SimulationBase<TSimType>::addHDF5Output()
   {
      BOOST_STATIC_ASSERT(sizeof(TSimType) == 0); 
   }

   template <typename TSimType> void SimulationBase<TSimType>::addASCIIOutput()
   {
      BOOST_STATIC_ASSERT(sizeof(TSimType) == 0); 
   }

   template <typename TSimType> void SimulationBase<TSimType>::initEquations()
   {
      BOOST_STATIC_ASSERT(sizeof(TSimType) == 0); 
   }

   template <typename TSimType> void SimulationBase<TSimType>::updateEquationsRTP(const int step)
   {
      BOOST_STATIC_ASSERT(sizeof(TSimType) == 0); 
   }

   template <typename TSimType> void SimulationBase<TSimType>::updateEquationsRHS()
   {
      BOOST_STATIC_ASSERT(sizeof(TSimType) == 0); 
   }

   template <typename TSimType> void SimulationBase<TSimType>::transformEquationsRHS(const int step)
   {
      BOOST_STATIC_ASSERT(sizeof(TSimType) == 0); 
   }

   template <typename TSimType> void SimulationBase<TSimType>::addExternalInfluence()
   {
      BOOST_STATIC_ASSERT(sizeof(TSimType) == 0); 
   }

   template <typename TSimType> void SimulationBase<TSimType>::timestepEquations()
   {
      BOOST_STATIC_ASSERT(sizeof(TSimType) == 0); 
   }

}

#endif // SIMULATIONBASE_HPP
