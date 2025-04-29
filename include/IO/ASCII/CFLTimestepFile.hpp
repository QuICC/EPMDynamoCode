/** \file CFLTimestepFile.hpp 
 *  \brief Implementation of a time information file
 */

#ifndef CFLTIMESTEPFILE_HPP
#define CFLTIMESTEPFILE_HPP

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
#include "IO/ASCII/ASCIIEWriter.hpp"
#include "Timestepping/TimestepParameters.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of ASCIIFile for a general timestep evolution ASCII file
    */
   class CFLTimestepFile: public ASCIIEWriter
   {
      public:
         /**
         * @brief Constructor
         *
         * @param name Name of the file
         * @param tsParams Timestep parameters
         */
         CFLTimestepFile(std::string name, const TimestepParameters &tsParams);

         /**
         * @brief Destructor
         */
         virtual ~CFLTimestepFile() {};

         /**
          * @brief Write spectrum to file
          */
         virtual void write();
         
      protected:

      private:

         /**
          * @brief Reference to a TimestepParameters object to get access to time and timestep
          */
         const TimestepParameters&  mrTSParams;
   };

   /// Typedef for a smart pointer of a CFLTimestepFile
   typedef EPMSHARED_PTR<CFLTimestepFile> SmartCFLTimestepFile;

}

#endif // CFLTIMESTEPFILE_HPP
