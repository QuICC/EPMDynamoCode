/** \file IOSystem.hpp
 *  \brief Implementation of the general control instance for the IO related operations
 */

#ifndef IOSYSTEM_HPP
#define IOSYSTEM_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "IO/IOSystemBase.hpp"
#include "IO/ASCII/ASCIIWriter.hpp"
#include "IO/HDF5/HDF5Writer.hpp"
#include "IO/HDF5/State/StateFileReaderBase.hpp"
#include "IO/HDF5/Imposed/ImposedFieldReaderBase.hpp"
#include "IO/HDF5/Source/SourceFileReaderBase.hpp"
#include "Timestepping/TimestepParameters.hpp"

namespace EPMDynamo {

   /**
    * @brief Full IO system reponsible for dealing with all the files
    */
   class IOSystem: public IOSystemBase
   {
      public:
         /**
         * @brief Constructor
         */
         IOSystem();

         /**
         * @brief Destructor
         */
         virtual ~IOSystem() {};

         /**
          * @brief Initialise the writer files
          */
         void initWriters();

         /**
          * @brief Finalise the writer files
          */
         void finaliseWriters();

         /**
          * @brief Set and use the given input statefile
          *
          * @param pInState Pointer to input state file
          * @param tsParams Timestep parameters
          */
         void useInitialState(SmartStateFileReaderBase pInState, TimestepParameters& tsParams);

         /**
          * @brief Use the given input imposefield file
          *
          * @param pIField Pointer to imposed field state file
          */
         void useImposedField(SmartImposedFieldReaderBase pIField);

         /**
          * @brief Use the given input source file
          *
          * @param pSource Pointer to source field state file
          */
         void useSource(SmartSourceReaderBase pSource);

         /**
          * @brief Add an output ASCIIWriter
          *
          * @param pWriter Pointer to ASCII writer
          */
         void addASCIIWriter(SmartASCIIWriter  pWriter);

         /**
          * @brief Add an output HDF5Writer
          *
          * @param pWriter Pointer to HDF5 writer
          */
         void addHDF5Writer(SmartHDF5Writer  pWriter);

         /**
          * @brief Write ASCII data
          */
         void writeASCII();

         /**
          * @brief Write HDF5 data
          */
         void writeHDF5();
         
      protected:

      private:
         /**
          * @brief Initialise the system
          */
         void init();

         /**
          * @brief Vector of ASCII output files
          */
         std::vector<SmartASCIIWriter> mASCIIWriters;

         /**
          * @brief Vector of HDF5 output files
          */
         std::vector<SmartHDF5Writer> mHDF5Writers;
   };

}

#endif // IOSYSTEM_HPP
