/** \file ControlInterface.cpp
 *  \brief Source of the external control interface implementation
 */

// Configuration includes
//
#include "Config/Parallelisation.h"

// System includes
//
#include <iostream>

// External includes
//

// Class include
//
#include "IO/ASCII/ControlInterface.hpp"

// Project includes
//
#include "IO/ASCII/ControlInterfaceDefs.hpp"
#include "General/EPMException.hpp"

namespace EPMDynamo {

   ControlInterface::ControlInterface()
      : ASCIIFile(ControlInterfaceDefs::NAME, ControlInterfaceDefs::EXTENSION, ControlInterfaceDefs::HEADER, ControlInterfaceDefs::VERSION), mKeepRunning(true)
   {
      this->init();
   }

   ControlInterface::~ControlInterface()
   {
      this->finalise();
   }

   void ControlInterface::init()
   {
      if(this->doesIO())
      {
         // Create file
         this->create();

         // Close file
         this->close();
      }
   }

   void ControlInterface::read()
   {
      int status = ControlInterfaceDefs::STATUS_GOON;

      if(this->doesIO())
      {
         // Open file
         status = this->open();

         // Close file
         this->close();
      }

      #ifdef EPMDYNAMO_MPI
         MPI_Bcast(&status, 1, MPI_INT, this->ioRank(), MPI_COMM_WORLD);
      #endif // EPMDYNAMO_MPI

      // Handle input and share the information
      this->analyseStatus(status);
   }

   void ControlInterface::finalise()
   {
      if(this->doesIO())
      {
         // Close the file
         this->close();

         // Delete the file
         this->deleteFile();
      }
   }

   void ControlInterface::create()
   {
      if(this->doesIO())
      {
         this->mFile.open(this->filename().c_str(), std::fstream::out);

         if(! this->mFile.is_open())
         {
            throw EPMException("ControlInterface::open", "Couldn't open file " + this->filename() + "!");
         }
      }
   }

   int ControlInterface::open()
   {
      this->mFile.open(filename().c_str());

      if(! this->mFile.is_open())
      {
         return ControlInterfaceDefs::STATUS_STOP;
      } else
      {
         return ControlInterfaceDefs::STATUS_GOON;
      }
   }

   void ControlInterface::close()
   {
      if(this->doesIO())
      {
         this->mFile.close();
      }
   }

   void ControlInterface::deleteFile()
   {
      std::remove(this->filename().c_str());
   }

   void ControlInterface::analyseStatus(const int status)
   {
      if(status == ControlInterfaceDefs::STATUS_STOP)
      {
         this->mKeepRunning = false;
         if(this->doesIO())
         {
            std::cout << "###########################" << std::endl;
            std::cout << "## User requested abort! ##" << std::endl;
            std::cout << "###########################" << std::endl;
         }
      }
   }

}
