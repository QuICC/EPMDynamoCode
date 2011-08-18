/** \file StateAnalyser.hpp
 *  \brief Implementation of a general state file analyser
 */

#ifndef STATEANALYSER_HPP
#define STATEANALYSER_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "IO/HDF5/State/StateFileReader.hpp"
#include "Utilities/GeneratorBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a general state file analyser
    *
    * \tparam TSimTraits Traits for the analyser
    */
   template <typename TSimTraits, template <typename> class TAnalysis> class StateAnalyser: public TAnalysis<TSimTraits>
   {
      public:
         /// Typedef for the state file reader
         typedef EPMSHARED_PTR<StateFileReader<TSimTraits> >  SmartStateReader; 

         /**
          * @brief Constructor
          */
         StateAnalyser();

         /**
          * @brief Simple empty destructor
          */
         virtual ~StateAnalyser() {};

         /**
          * @brief Initialise the output file
          *
          * @param name Basename of the file
          */
         void setup(std::string name = "2Analyse");

         /**
          * @brief Finalise the output cscs file
          */
         void finalise();

      protected:

         /**
          * @brief Smart pointer for the state file reader
          */
         SmartStateReader  mpInFile;

         /**
          * @brief Initialise the input file
          */
         void initInputFile(std::string name);

      private:
   };

   template <typename TSimTraits, template <typename> class TAnalysis> StateAnalyser<TSimTraits,TAnalysis>::StateAnalyser()
      : TAnalysis<TSimTraits>()
   {
   }

   template <typename TSimTraits, template <typename> class TAnalysis> void StateAnalyser<TSimTraits,TAnalysis>::setup(std::string inName)
   {
      // Init the input file
      this->initInputFile(inName);
   }

   template <typename TSimTraits, template <typename> class TAnalysis> void StateAnalyser<TSimTraits,TAnalysis>::initInputFile(std::string inName)
   {
      SmartStateReader  pInFile;

      // Create file with all fields
      if(TSimTraits::NeedCodensity && TSimTraits::NeedMagnetic && TSimTraits::NeedVelocity)
      {
         pInFile.reset(new StateFileReader<TSimTraits>(this->codC(), this->magB(), this->velV(), inName));

      // Create file with Codensity and Velocity fields
      } else if(TSimTraits::NeedCodensity && TSimTraits::NeedVelocity)
      {
         pInFile.reset(new StateFileReader<TSimTraits>(this->codC(), this->velV(), inName));

      // Create file with Magnetic and Velocity fields
      } else if(TSimTraits::NeedMagnetic && TSimTraits::NeedVelocity)
      {
         pInFile.reset(new StateFileReader<TSimTraits>(this->magB(), this->velV(), inName));

      // Create file with only codensity field
      } else if(TSimTraits::NeedCodensity)
      {
         pInFile.reset(new StateFileReader<TSimTraits>(this->codC(), inName));

      // Create file with only magnetic field
      } else if(TSimTraits::NeedMagnetic)
      {
         pInFile.reset(new StateFileReader<TSimTraits>(this->magB(), inName));

      // Create file with only velocity field
      } else if(TSimTraits::NeedVelocity)
      {
         pInFile.reset(new StateFileReader<TSimTraits>(this->velV(), inName));
      }
      // Set the state file localy
      this->mpInFile = pInFile;

      // Initialise the input file
      this->mpInFile->init();

      // Read setup information from file
      this->mpInFile->read();
   }

   template <typename TSimTraits, template <typename> class TAnalysis> void StateAnalyser<TSimTraits,TAnalysis>::finalise()
   {
      // Finalise input file IO
      this->mpInFile->finalise();

      // Finalise generator base
      GeneratorBase<TSimTraits>::finalise();
   }

}

#endif // STATEANALYSER_HPP
