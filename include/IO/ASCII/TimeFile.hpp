/** \file TimeFile.hpp 
 *  \brief Implementation of a time information file
 */

#ifndef TIMEFILE_HPP
#define TIMEFILE_HPP

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
   class TimeFile: public ASCIIEWriter
   {
      public:
         /**
         * @brief Constructor
         *
         * @param name Name of the file
         * @param tsParams Timestep parameters
         */
         TimeFile(std::string name, const TimestepParameters &tsParams);

         /**
         * @brief Destructor
         */
         virtual ~TimeFile() {};

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

   /// Typedef for a smart pointer of a TimeFile
   typedef EPMSHARED_PTR<TimeFile> SmartTimeFile;

}

#endif // TIMEFILE_HPP
