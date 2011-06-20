/** \file CodensityDiffusionSimulation.hpp
 *  \brief Implementation of a codensity diffusion simulation
 */

#ifndef CODENSITYDIFFUSIONSIMULATION_HPP
#define CODENSITYDIFFUSIONSIMULATION_HPP

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
#include "Simulations/SimulationBase.hpp"
#include "Simulations/Traits/CodensityDiffusionTraits.hpp"

#include "IO/HDF5/State/StateFileReader.hpp"
#include "IO/HDF5/State/StateFileWriter.hpp"
#include "IO/HDF5/Source/SourceFileReader.hpp"
#include "IO/ASCII/EnergyFile.hpp"
#include "IO/ASCII/SpectrumFile.hpp"
#include "IO/ASCII/TimeFile.hpp"

#include "Equations/Transport/TransportDiffusion.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a codensity diffusion simulation
    */
   class CodensityDiffusionSimulation: public SimulationBase
   {
      public:
         /**
          * @brief Simple empty destructor
          */
         virtual ~CodensityDiffusionSimulation() {};

      protected:
         /**
          * @brief Constructor
          */
         CodensityDiffusionSimulation();

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
          * \param step Current step in a multistep transform
          */
         void updateEquationsRTP(const int step);

         /**
          * @brief Update equations RHS
          */
         void updateEquationsRHS();

         /**
          * @brief Transform equations RHS
          *
          * \param step Current step in a multistep transform
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

         /**
          * @brief Configure the transforms
          */
         void configureTransforms();

      private:
         /**
          * @brief Codensity scalar
          */
         CodensityDiffusionTraits::CodType   mCodC;

         /**
          * @brief Transport equation
          */
         TransportDiffusion<CodensityDiffusionTraits>    mTransport;
   };

}

#endif // CODENSITYDIFFUSIONSIMULATION_HPP
