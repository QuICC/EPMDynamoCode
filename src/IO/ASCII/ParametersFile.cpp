/** \file ParametersFile.cpp
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
#include "IO/ASCII/ParametersFile.hpp"

// Project includes
//
#include "General/EPMException.hpp"
#include "IO/ASCII/ParametersFileDefs.hpp"

namespace EPMDynamo {

   ParametersFile::ParametersFile(std::string type)
      : XMLReader(type, ParametersFileDefs::BASENAME, ParametersFileDefs::EXTENSION, ParametersFileDefs::HEADER, ParametersFileDefs::VERSION), mTruncArray(5), mEqArray(1), mTStepArray(2), mRunArrayI(3), mRunArray(1)
   {
   }

   void ParametersFile::read()
   {
      try{
         if(this->doesIO())
         {
            // Read truncation values from file
            this->readTruncation();

            // Read Physical parameters from file
            this->readPhysical();

            // Read timestepping parameters from file
            this->readTimestepping();

            // Read run parameters from file
            this->readRun();
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

   std::string ParametersFile::stdMessageName() const
   {
      std::string fname = "OUT";
      std::ostringstream  oss;

      oss << this->mTruncArray(0);
      fname += "N" + oss.str();
      oss.str("");

      oss << this->mTruncArray(1);
      fname += "L" + oss.str();
      oss.str("");

      oss << this->mTruncArray(2);
      fname += "M" + oss.str();
      oss.str("");

      oss << this->mTruncArray(3);
      fname += "Mp" + oss.str();
      oss.str("");

      return fname;
   }

   const ArrayI& ParametersFile::aTrunc() const
   {
      return this->mTruncArray;
   }

   const Array& ParametersFile::aTStep() const
   {
      return this->mTStepArray;
   }

   const ArrayI& ParametersFile::aRunI() const
   {
      return this->mRunArrayI;
   }
 
   const Array& ParametersFile::aRun() const
   {
      return this->mRunArray;
   }

   const Array& ParametersFile::aEquation() const
   {
      return this->mEqArray;
   }

   void ParametersFile::readTruncation()
   {
      if(this->doesIO())
      {
         rapidxml::xml_node<> *node = mXML.first_node(ParametersFileDefs::TRUNCXML.c_str());

         if(node)
         {
            this->readValue(this->mTruncArray(0), node, ParametersFileDefs::TRUNCNXML);

            this->readValue(this->mTruncArray(1), node, ParametersFileDefs::TRUNCLXML);

            this->readValue(this->mTruncArray(2), node, ParametersFileDefs::TRUNCMXML);

            this->readValue(this->mTruncArray(3), node, ParametersFileDefs::TRUNCMPXML);

            this->readValue(this->mTruncArray(4), node, ParametersFileDefs::TRUNCNCOREXML);
         } else
         {
            throw EPMException("ParametersFile::readTruncation", "Couldn't find tag!");
         }

         if(this->mTruncArray(1) < this->mTruncArray(2))
         {
            throw EPMException("ParametersFile::readTruncation", "Bad truncation: M > L");
         }
      } else
      {
         this->mTruncArray.setConstant(-1);
      }
   }

   void ParametersFile::readPhysical()
   {
      if(this->doesIO())
      {
         rapidxml::xml_node<> *node = mXML.first_node(ParametersFileDefs::PHYSICALXML.c_str());

         if(node)
         {
            if(this->mType == "EkQRaRo")
            {
               // Resize equation array to correct size
               this->mEqArray.resize(4);

               // Read the Ekman number
               this->readValue(this->mEqArray(0), node, ParametersFileDefs::PHYSEKMANXML);

               // Read the Roberts number
               this->readValue(this->mEqArray(1), node, ParametersFileDefs::PHYSROBERTSXML);

               // Read the Rayleigh number
               this->readValue(this->mEqArray(2), node, ParametersFileDefs::PHYSRAYLEIGHXML);

               // Read the Rossby number
               this->readValue(this->mEqArray(3), node, ParametersFileDefs::PHYSROSSBYXML);
            } else
            {
               throw EPMException("ParametersFile::readPhysical", "The requested type is not implemented! (yet?)");
            }
         } else
         {
            throw EPMException("ParametersFile::readPhysical", "Couldn't find tag!");
         }
      } else
      {
         this->mEqArray.setConstant(-1);
      }
   }

   void ParametersFile::readTimestepping()
   {
      if(this->doesIO())
      {
         rapidxml::xml_node<> *node = mXML.first_node(ParametersFileDefs::TIMEXML.c_str());

         if(node)
         {
            this->readValue(this->mTStepArray(0), node, ParametersFileDefs::TIMETIMEXML);

            this->readValue(this->mTStepArray(1), node, ParametersFileDefs::TIMETSTEPXML);
         } else
         {
            throw EPMException("ParametersFile::readTimestepping", "Couldn't find tag!");
         }
      } else
      {
         this->mTStepArray.setConstant(-1);
      }
   }

   void ParametersFile::readRun()
   {
      if(this->doesIO())
      {
         rapidxml::xml_node<> *node = mXML.first_node(ParametersFileDefs::RUNXML.c_str());

         if(node)
         {
            this->readValue(this->mRunArrayI(0), node, ParametersFileDefs::RUNMAXTSTEPXML);

            this->readValue(this->mRunArrayI(1), node, ParametersFileDefs::RUNARATEXML);

            this->readValue(this->mRunArrayI(2), node, ParametersFileDefs::RUNSRATEXML);

            this->readValue(this->mRunArray(0), node, ParametersFileDefs::RUNWALLXML);
         } else
         {
            throw EPMException("ParametersFile::readRun", "Couldn't find tag!");
         }
      } else
      {
         this->mRunArray.setConstant(-1);
         this->mRunArrayI.setConstant(-1);
      }
   }

   void ParametersFile::printInfo() const
   {
      if(this->doesIO())
      {
         std::cout << "--------------------" << std::endl;
         std::cout << "**** Truncation ****" << std::endl;
         std::cout << "--------------------" << std::endl;
         std::cout << "  " << "N: " << this->mTruncArray(0) << std::endl;
         std::cout << "  " << "L: " << this->mTruncArray(1) << std::endl;
         std::cout << "  " << "M: " << this->mTruncArray(2) << std::endl;
         std::cout << "  " << "Mp: " << this->mTruncArray(3) << std::endl;
         std::cout << "  " << "Ncore: " << this->mTruncArray(4) << std::endl;
         std::cout << "********************" << std::endl;
         std::cout << std::endl;
         std::cout << "--------------------" << std::endl;
         std::cout << "**** Equations *****" << std::endl;
         std::cout << "--------------------" << std::endl;
         if(this->mType == "EkQRaRo")
         {
            std::cout << "  " << "E: " << this->mEqArray(0) << std::endl;
            std::cout << "  " << "q: " << this->mEqArray(1) << std::endl;
            std::cout << "  " << "Ra: " << this->mEqArray(2) << std::endl;
            std::cout << "  " << "Ro: " << this->mEqArray(3) << std::endl;
         } else
         {
            std::cout << " !!!! UNKNOWN PHYSICAL PARAMETERS !!!! " << std::endl;
         }
         std::cout << "********************" << std::endl;
         std::cout << std::endl;
         std::cout << "--------------------" << std::endl;
         std::cout << "*** Timestepping ***" << std::endl;
         std::cout << "--------------------" << std::endl;
         std::cout << "  " << "Time: " << this->mTStepArray(0) << std::endl;
         std::cout << "  " << "Timestep: " << this->mTStepArray(1) << std::endl;
         std::cout << "********************" << std::endl;
         std::cout << std::endl;
         std::cout << "--------------------" << std::endl;
         std::cout << "******* Run ********" << std::endl;
         std::cout << "--------------------" << std::endl;
         std::cout << "  " << "Maxtstep: " << this->mRunArrayI(0) << std::endl;
         std::cout << "  " << "ASCII Rate: " << this->mRunArrayI(1) << std::endl;
         std::cout << "  " << "State Rate: " << this->mRunArrayI(2) << std::endl;
         std::cout << "  " << "Walltime: " << this->mRunArray(0) << std::endl;
         std::cout << "********************" << std::endl;
         std::cout << std::endl;
      }
   }

   void ParametersFile::spreadParameters()
   {
      #ifdef EPMDYNAMO_MPI

      int nBlocks = 5;
      int idx = 0;
      MPI_Aint    displ[nBlocks];
      int         blocks[nBlocks];
      MPI_Datatype   types[nBlocks];

      MPI_Aint    element;

      // Create mTruncArray part
      MPI_Get_address(this->mTruncArray.data(), &element);
      displ[idx] = element;
      blocks[idx] = 5;
      types[idx] = MPI_INT;
      idx++;

      // Create mEqArray part
      MPI_Get_address(this->mEqArray.data(), &element);
      displ[idx] = element;
      blocks[idx] = this->mEqArray.size();
      types[idx] = MPI_DOUBLE;
      idx++;

      // Create mTStepArray part
      MPI_Get_address(this->mTStepArray.data(), &element);
      displ[idx] = element;
      blocks[idx] = 2;
      types[idx] = MPI_DOUBLE;
      idx++;

      // Create mRunArrayI part
      MPI_Get_address(this->mRunArrayI.data(), &element);
      displ[idx] = element;
      blocks[idx] = 3;
      types[idx] = MPI_INT;
      idx++;

      // Create mRunArray part
      MPI_Get_address(this->mRunArray.data(), &element);
      displ[idx] = element;
      blocks[idx] = 1;
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
