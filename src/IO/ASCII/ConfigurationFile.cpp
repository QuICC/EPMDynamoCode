/** \file ConfigurationFile.cpp
 *  \brief Source of the implementation of the XML parameters file reader
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
#include "IO/ASCII/ConfigurationFile.hpp"

// Project includes
//
#include "General/EPMException.hpp"
#include "IO/ASCII/ConfigurationFileDefs.hpp"
#include "Timestepping/TimestepConfig.hpp"

namespace EPMDynamo {

   ConfigurationFile::ConfigurationFile(const std::string type, const std::vector<std::string>& integers, const std::vector<std::string>& floats)
      : XMLReader(type, ConfigurationFileDefs::BASENAME, ConfigurationFileDefs::EXTENSION, ConfigurationFileDefs::HEADER, ConfigurationFileDefs::VERSION), mIntegerIDs(integers), mFloatIDs(floats)
   {
      // Resize arrays depending
      this->setupStorage();
   }

   void ConfigurationFile::setupStorage()
   {
      // Resize integer array
      this->mIntegers.resize(this->mIntegerIDs.size());

      // Resize float array
      this->mFloats.resize(this->mFloatIDs.size());
   }

   void ConfigurationFile::read()
   {
      try{
         if(this->doesIO())
         {
            // Read integer values from file
            this->readIntegers();

            // Read floating values from file
            this->readFloats();
         }
         
         // Share Parameters with other CPUs (if applicatble)
         this->spreadParameters();
      }
      catch(EPMException &e)
      {
         e.printStdMessage();

         throw -1;
      }
   }

   void ConfigurationFile::readIntegers()
   {
      if(this->doesIO())
      {
         rapidxml::xml_node<> *node = this->mXML.first_node(ConfigurationFileDefs::CONFIGXML.c_str());

         if(node)
         {
            for(int i = 0; i < this->mIntegers.size(); i++)
            {
               this->readValue(this->mIntegers(i), node, this->mIntegerIDs.at(i));
            }
         } else
         {
            throw EPMException("ConfigurationFile::readIntegers", "Couldn't find tag!");
         }
      } else
      {
         this->mIntegers.setConstant(-1);
      }
   }

   void ConfigurationFile::readFloats()
   {
      if(this->doesIO())
      {
         rapidxml::xml_node<> *node = mXML.first_node(ConfigurationFileDefs::CONFIGXML.c_str());

         if(node)
         {
            for(int i = 0; i < this->mFloats.size(); i++)
            {
               this->readValue(this->mFloats(i), node, this->mFloatIDs.at(i));
            }
         } else
         {
            throw EPMException("ConfigurationFile::readFloats", "Couldn't find tag!");
         }
      } else
      {
         this->mFloats.setConstant(-1);
      }
   }

   void ConfigurationFile::printInfo() const
   {
      if(this->doesIO())
      {
         std::cout << "--------------------" << std::endl;
         std::cout << "**** Configuration ****" << std::endl;
         std::cout << "--------------------" << std::endl;
         for(int i = 0; i < this->mIntegers.size(); i++)
         {
            std::cout << "  " << this->mIntegerIDs.at(i) << ": " << this->mIntegers(i) << std::endl;
         }
         for(int i = 0; i < this->mFloats.size(); i++)
         {
            std::cout << "  " << this->mFloatIDs.at(i) << ": " << this->mFloats(i) << std::endl;
         }
         std::cout << "********************" << std::endl;
         std::cout << std::endl;
      }
   }

   void ConfigurationFile::spreadParameters()
   {
      #ifdef EPMDYNAMO_MPI

      int nBlocks = 2;
      int idx = 0;
      MPI_Aint    displ[nBlocks];
      int         blocks[nBlocks];
      MPI_Datatype   types[nBlocks];

      MPI_Aint    element;

      // Create integer part
      MPI_Get_address(this->mIntegers.data(), &element);
      displ[idx] = element;
      blocks[idx] = this->mIntegers.size();
      types[idx] = MPI_INT;
      idx++;

      // Create floats part
      MPI_Get_address(this->mFloats.data(), &element);
      displ[idx] = element;
      blocks[idx] = this->mFloats.size();
      types[idx] = MPI_DOUBLE;
      idx++;

      MPI_Datatype   paramType;

      // Create MPI datatype
      MPI_Type_create_struct(nBlocks, blocks, displ, types, &paramType);
      // Commit MPI datatype
      MPI_Type_commit(&paramType);

      // Broadcast the information
      MPI_Bcast(MPI_BOTTOM, 1, paramType, this->ioRank(), MPI_COMM_WORLD);
      
      // Free the datatype
      MPI_Type_free(&paramType);

      #endif // EPMDYNAMO_MPI
   }
}
