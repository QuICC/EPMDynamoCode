/** \file Simulation.hpp
 *  \brief High level implementation of a simulation
 */

#ifndef SIMULATION_HPP
#define SIMULATION_HPP

// Configuration includes
//
#include "Config/Parallelisation.h"

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "General/EPMException.hpp"

namespace EPMDynamo {

   /**
    * @brief High level implementation of a simulation's execution steps.
    *
    * The setup and actual definition of the different steps are done in the TSimImpl object.
    *
    * \tparam TSimType Type of simulation
    * \tparam TSimImpl Actual implementation of the simulation
    */
   template <typename TSimType, template<typename > class TSimImpl> class Simulation: public TSimImpl<TSimType>
   {
      public:
         /**
          * @brief Constructor
          *
          * The constructor simply calls the constructor of TSimImpl.
          */
         Simulation();

         /**
          * @brief Simple empty destructor
          */
         virtual ~Simulation() {};

         /**
          * @brief Initialise the different components of the simulation
          */
         void init();

         /**
          * @brief Run the simulation
          */
         void run();

         /**
          * @brief Finalise simulation run
          */
         void finalise();

      protected:

      private:
   };

   template <typename TSimType, template<typename > class TSimImpl> Simulation<TSimType, TSimImpl>::Simulation()
      : TSimImpl<TSimType>()
   {
      // Stop timer for construction time
      this->mExecTimer.stop();
      this->mExecTimer.update(1);

      // Print some information on simulation construction time
      if(this->mpTrunc->para().id() == 0)
      {
         std::cout << std::endl;
         std::cout << "-------------------------------------------" << std::endl;
         std::cout << "*** Simulation construction successfull ***" << std::endl;
         std::cout << "-------------------------------------------" << std::endl;
         std::cout << std::endl;
      }
   }

   template <typename TSimType, template<typename > class TSimImpl> void Simulation<TSimType, TSimImpl>::init()
   {
      // Start timer
      this->mExecTimer.start();

      try{
         // Initialise the simulation equations
         this->initEquations();

         // Configure the real space <-> spectral space transforms
         this->configureTransforms();

         // Initialise the simulation fields
         this->initFields();

         // Add the Ascii diagnostics output
         this->addASCIIOutput();

         // Add the HDF5 output
         this->addHDF5Output();

         // Initialise the ouputs
         this->initOutput();
      }
      catch(EPMException &e)
      {
         e.printStdMessage();

         throw -1;
      }

      // Make sure CPUs are synchronized after initialisation
      EPMDYNAMO_SYNCHRONIZE;

      // Stop timer
      this->mExecTimer.stop();
      this->mExecTimer.update(2);

      // Produce some nice output information for initialisation step
      if(this->mpTrunc->para().id() == 0)
      {
         std::cout << "---------------------------------------------" << std::endl;
         std::cout << "*** Simulation initialisation successfull ***" << std::endl;
         std::cout << "---------------------------------------------" << std::endl;
         std::cout << std::endl;
      }
   }

   template <typename TSimType, template<typename > class TSimImpl> void Simulation<TSimType, TSimImpl>::run()
   {
      // Start timer
      this->mExecTimer.start();

      // Execute last initialisation before run
      this->preRun();

      // Stop pre run timing
      this->mExecTimer.stop();
      this->mExecTimer.update(3);

      // Start timer
      this->mExecTimer.start();

      // Start looping as long as requested
      while(this->mSimControl.keepRunning())
      {
         // Update RHS of all equation
         for(int i=0; i < this->mTransformSteps; ++i)
         {
            this->combineRTPTransforms(i);

            this->updateEquationsRTP(i);
         }

         // Update RHS of all equation
         this->updateEquationsRHS();

         // Transform RHS to spectral space
         for(int i=0; i < this->mTransformSteps; ++i)
         {
            this->combineSpectralTransforms(i);

            this->transformEquationsRHS(i);
         }

         // Include possible source terms / imposed fields
         this->addExternalInfluence();

         // Update the timestep
         this->mSimControl.tsControl().updateTimestep();

         // Timestep all equations
         this->timestepEquations();

         // Update other values of simulation control
         this->mSimControl.update();

         // Save files
         this->writeFiles();

         // Synchronize simulation
         EPMDYNAMO_SYNCHRONIZE;
      }
      // Stop pre run timing
      this->mExecTimer.stop();
      this->mExecTimer.update(4);

      // Start timer
      this->mExecTimer.start();

      // Execute last operations before end
      this->postRun();

      // Synchronize simulation
      EPMDYNAMO_SYNCHRONIZE;

      // Stop post run timing
      this->mExecTimer.stop();
      this->mExecTimer.update(5);
   }

   template <typename TSimType, template<typename > class TSimImpl> void Simulation<TSimType, TSimImpl>::finalise()
   {
      // Print timestepping infos
      this->mSimControl.printInfo();

      // Print execution timer infos
      this->mExecTimer.analyze();
      this->mExecTimer.printInfo();

      // Close and finalise writers
      this->mIOSys.finaliseWriters();
   }

}

#endif // SIMULATION_HPP
