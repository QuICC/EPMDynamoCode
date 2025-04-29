/** \file ParametersFile.hpp 
 *  \brief Implementation of the parameters file reader
 */

#ifndef PARAMETERSFILE_HPP
#define PARAMETERSFILE_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "IO/ASCII/XMLReader.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of the XML parameters file reader
    */
   class ParametersFile: public XMLReader
   {
      public:
         /**
         * @brief Constructor
         *
         * @param type Type of the parameters file
         */
         ParametersFile(std::string type);

         /**
         * @brief Destructor
         */
         virtual ~ParametersFile() {};

         /**
          * @brief Get truncation parameters array
          */
         const ArrayI& aTrunc() const;

         /**
          * @brief Get timestep parameters array
          */
         const Array& aTStep() const;

         /**
          * @brief Get integer run parameters array
          */
         const ArrayI& aRunI() const;

         /**
          * @brief Get run parameters array
          */
         const Array& aRun() const;

         /**
          * @brief Get equation parameters array
          */
         const Array& aEquation() const;

         /**
          * @brief Get boundary condition parameters array
          */
         const ArrayI& aBC() const;

         /**
          * @brief Read content of configuration file
          */
         virtual void read();

         /**
          * @brief Create StdMessage file name
          */
         std::string stdMessageName() const;

         /**
          * @brief Output run information
          */
         void printInfo() const;
         
      protected:
         /**
          * @brief Resize arrays depending on file type
          */
         void setupStorage();

         /**
          * @brief Read truncation information from file
          */
         void readTruncation();

         /**
          * @brief Read physical parameters information from file
          */
         void readPhysical();

         /**
          * @brief Read boundary condition parameters information from file
          */
         void readBoundary();

         /**
          * @brief Read io parameters information from file
          */
         void readTimestepping();

         /**
          * @brief Read run parameters information from file
          */
         void readRun();

      private:
         /**
          * @brief Storage for the truncation values
          */
         ArrayI   mTruncArray;

         /**
          * @brief Storage for the equation parameters
          */
         Array    mEqArray;

         /**
          * @brief Storage for the boundary paramters
          */
         ArrayI    mBCArray;

         /**
          * @brief Storage for the timestep parameters
          */
         Array    mTStepArray;

         /**
          * @brief Storage for the integer run parameters
          */
         ArrayI    mRunArrayI;

         /**
          * @brief Storage for the run parameters
          */
         Array     mRunArray;

         /**
          * @brief Spread parameters over parallel simulation
          */
         void spreadParameters();
   };

   /// Typedef for a smart pointer of a Parameters
   typedef EPMSHARED_PTR<ParametersFile> SmartParametersFile;

}

#endif // PARAMETERSFILE_HPP
