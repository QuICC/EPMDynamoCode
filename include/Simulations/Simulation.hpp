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
    * The setup and actual definition of the different steps are done in by the TSimImpl object.
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
      this->mTimer.stop();

      // Print some information on simulation construction time
      if(this->mpTrunc->para().id() == 0)
      {
         std::cout << std::endl;
         std::cout << "----------------------------------------------------------" << std::endl;
         std::cout << "*** Simulation construction successfull in " << static_cast<int>(std::ceil(this->mTimer.time())) << " seconds ***" << std::endl;
         std::cout << "----------------------------------------------------------" << std::endl;
         std::cout << std::endl;
      }
   }

   template <typename TSimType, template<typename > class TSimImpl> void Simulation<TSimType, TSimImpl>::init()
   {
      // Start timer
      this->mTimer.start();

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
      this->mTimer.stop();

      // Produce some nice output information for initialisation step
      if(this->mpTrunc->para().id() == 0)
      {
         std::cout << "------------------------------------------------" << std::endl;
         if(this->mTimer.time() < 1.0)
         {
            std::cout << "*** Initialisation successfull in < 1 second ***" << std::endl;
         } else
         {
            std::cout << "*** Initialisation successfull in " << static_cast<int>(std::ceil(this->mTimer.time())) << " seconds ***" << std::endl;
         }
         std::cout << "------------------------------------------------" << std::endl;
         std::cout << std::endl;
      }
   }

   template <typename TSimType, template<typename > class TSimImpl> void Simulation<TSimType, TSimImpl>::run()
   {
      // Execute last initialisation before run
      this->preRun();

EPMFloat t_upRTP(0.0);
EPMFloat t_upRHS(0.0);
EPMFloat t_trRHS(0.0);
EPMFloat t_exInf(0.0);
EPMFloat t_upTim(0.0);
EPMFloat t_tiEqu(0.0);
EPMFloat t_siUpd(0.0);
EPMFloat t_wrFil(0.0);
EPMFloat t_sync(0.0);

      // Start looping as long as requested
      while(this->mSimControl.keepRunning())
      {
this->mTimer.start();
         // Update RHS of all equation
         for(int i=0; i < this->mTransformSteps; ++i)
         {
            this->combineRTPTransforms(i);

            this->updateEquationsRTP(i);
         }
this->mTimer.stop();
t_upRTP += this->mTimer.time();

this->mTimer.start();
         // Update RHS of all equation
         this->updateEquationsRHS();
this->mTimer.stop();
t_upRHS += this->mTimer.time();

this->mTimer.start();
         // Transform RHS to spectral space
         for(int i=0; i < this->mTransformSteps; ++i)
         {
            this->combineSpectralTransforms(i);

            this->transformEquationsRHS(i);
         }
this->mTimer.stop();
t_trRHS += this->mTimer.time();

this->mTimer.start();
         // Include possible source terms / imposed fields
         this->addExternalInfluence();
this->mTimer.stop();
t_exInf += this->mTimer.time();

this->mTimer.start();
         // Update the timestep
         this->mSimControl.tsControl().updateTimestep();
this->mTimer.stop();
t_upTim += this->mTimer.time();

this->mTimer.start();
         // Timestep all equations
         this->timestepEquations();
this->mTimer.stop();
t_tiEqu += this->mTimer.time();

this->mTimer.start();
         // Update other values of simulation control
         this->mSimControl.update();
this->mTimer.stop();
t_siUpd += this->mTimer.time();

this->mTimer.start();
         // Save files
         this->writeFiles();
this->mTimer.stop();
t_wrFil += this->mTimer.time();

this->mTimer.start();
         // Synchronize simulation
         EPMDYNAMO_SYNCHRONIZE;
this->mTimer.stop();
t_sync += this->mTimer.time();
      }

std::ofstream  file;
std::ostringstream oss;
oss << EPMDYNAMO_RANK;
file.open("sim_timings"+oss.str()+".dat");
file << t_upRTP << std::endl;
file << t_upRHS << std::endl;
file << t_trRHS << std::endl;
file << t_exInf << std::endl;
file << t_upTim << std::endl;
file << t_tiEqu << std::endl;
file << t_siUpd << std::endl;
file << t_wrFil << std::endl;
file << t_sync << std::endl;

file << std::endl << std::endl << std::endl;
file << t_upRTP << "\t";
file << t_upRHS << "\t";
file << t_trRHS << "\t";
file << t_exInf << "\t";
file << t_upTim << "\t";
file << t_tiEqu << "\t";
file << t_siUpd << "\t";
file << t_wrFil << "\t";
file << t_sync << std::endl;
file.close();

      // Execute last initialisation before end
      this->postRun();

      // Synchronize simulation
      EPMDYNAMO_SYNCHRONIZE;
   }

   template <typename TSimType, template<typename > class TSimImpl> void Simulation<TSimType, TSimImpl>::finalise()
   {
      // Print timestepping infos
      this->mSimControl.printInfo();

      // Close and finalise writers
      this->mIOSys.finaliseWriters();
   }

}

#endif // SIMULATION_HPP
