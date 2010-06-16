/** \file ControlInterface.hpp
 *  \brief Implementation of a external runtime control interface
 */

#ifndef CONTROLINTERFACE_HPP
#define CONTROLINTERFACE_HPP

// System includes
//
#include <fstream>

// External includes
//

// Project includes
//
#include "IO/ASCII/ASCIIFile.hpp"

namespace EPMDynamo {

   /**
    * @brief Creates an ascii file for runtime user interaction (abort etc...)
    */
   class ControlInterface: public ASCIIFile
   {
      public:
         /**
         * @brief Constructor
         */
         ControlInterface();

         /**
         * @brief Destructor
         */
         virtual ~ControlInterface();

         /**
          * @brief Create control file file
          */
         virtual void init();

         /**
          * @brief Read input or check existance
          */
         virtual void read();

         /**
          * @brief Finalise the file
          */
         virtual void finalise();

         /**
          * @brief User input leaves simulation running?
          */
         bool keepRunning() const;
         
      protected:
         /**
          * @brief Handle to the file
          */
         std::fstream mFile;

         /**
          * @brief Create the file
          */
         void create();

         /**
          * @brief Open the file
          */
         int open();

         /**
          * @brief Close the file
          */
         void close();

         /**
          * @brief Delete the file
          */
         void deleteFile();

         /**
          * @brief Analyse input interface status
          */
         void analyseStatus(const int status);

      private:
         /**
          * @brief Keep simulation running?
          */
         bool mKeepRunning;
   };

   inline bool ControlInterface::keepRunning() const
   {
      return this->mKeepRunning;
   }
}

#endif // CONTROLINTERFACE_HPP
