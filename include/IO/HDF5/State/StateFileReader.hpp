/** \file StateFileReader.hpp
 *  \brief Implementation of HDF5 state file reader
 */

#ifndef STATEFILEREADER_HPP
#define STATEFILEREADER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "IO/HDF5/StateFileReaderBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of HDF5 state file reader
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class StateFileReader: public StateFileReaderBase
   {
      public:
         /**
         * @brief Constructor with all three fields
         *
         * @param codC Codensity variable
         * @param magB Magnetic variable
         * @param velV Velocity variable
         * @param name file name
         */
         StateFileReader(typename TSimTraits<TSimType>::CodType& codC, typename TSimTraits<TSimType>::MagType& magB, typename TSimTraits<TSimType>::VelType& velV, std::string name);

         /**
         * @brief Constructor for thermal convection
         *
         * @param codC Codensity variable
         * @param velV Velocity variable
         * @param name file name
         */
         StateFileReader(typename TSimTraits<TSimType>::CodType& codC, typename TSimTraits<TSimType>::VelType& velV, std::string name);

         /**
         * @brief Constructor for magneto-convection
         *
         * @param magB Magnetic variable
         * @param velV Velocity variable
         * @param name file name
         */
         StateFileReader(typename TSimTraits<TSimType>::MagType& magB, typename TSimTraits<TSimType>::VelType& velV, std::string name);

         /**
         * @brief Constructor with only codensity
         *
         * @param codC Codensity variable
         * @param name file name
         */
         StateFileReader(typename TSimTraits<TSimType>::CodType& codC, std::string name);

         /**
         * @brief Constructor with only magnetic field
         *
         * @param magB Magnetic variable
         * @param name file name
         */
         StateFileReader(typename TSimTraits<TSimType>::MagType& magB, std::string name);

         /**
         * @brief Constructor for cases with all three fields
         *
         * @param velV Velocity variable
         * @param name file name
         */
         StateFileReader(typename TSimTraits<TSimType>::VelType& velV, std::string name);

         /**
         * @brief Destructor
         */
         virtual ~StateFileReader() {};

         /**
          * @brief Read State to file
          */
         virtual void read();
         
      protected:
         /**
          * @brief Pointer to the condensity variable
          */
         typename TSimTraits<TSimType>::CodType*  mpCodC;

         /**
          * @brief Pointer to the magnetic variable
          */
         typename TSimTraits<TSimType>::MagType*  mpMagB;

         /**
          * @brief Pointer to the velocity variable
          */
         typename TSimTraits<TSimType>::VelType*  mpVelV;

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits> StateFileReader<TSimType, TSimTraits>::StateFileReader(typename TSimTraits<TSimType>::CodType &codC, typename TSimTraits<TSimType>::MagType &magB, typename TSimTraits<TSimType>::VelType &velV, std::string name)
      : StateFileReaderBase(name, codC.trunc()), mpCodC(codC), mpMagB(magB), mpVelV(velV)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> StateFileReader<TSimType, TSimTraits>::StateFileReader(typename TSimTraits<TSimType>::CodType &codC, typename TSimTraits<TSimType>::VelType &velV, std::string name)
      : StateFileReaderBase(name, codC.trunc()), mpCodC(codC), mpMagB(NULL), mpVelV(velV)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> StateFileReader<TSimType, TSimTraits>::StateFileReader(typename TSimTraits<TSimType>::MagType &magB, typename TSimTraits<TSimType>::VelType &velV, std::string name)
      : StateFileReaderBase(name, codC.trunc()), mpCodC(NULL), mpMagB(magB), mpVelV(velV)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> StateFileReader<TSimType, TSimTraits>::StateFileReader(typename TSimTraits<TSimType>::CodType &codC, std::string name)
      : StateFileReaderBase(name, codC.trunc()), mpCodC(codC), mpMagB(NULL), mpVelV(NULL)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> StateFileReader<TSimType, TSimTraits>::StateFileReader(typename TSimTraits<TSimType>::MagType &magB, std::string name)
      : StateFileReaderBase(name, magB.trunc()), mpCodC(NULL), mpMagB(magB), mpVelV(NULL)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> StateFileReader<TSimType, TSimTraits>::StateFileReader(typename TSimTraits<TSimType>::VelType &codC, std::string name)
      : StateFileReaderBase(name, velV.trunc()), mpCodC(NULL), mpMagB(NULL), mpVelV(velV)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void StateFileReader<TSimType, TSimTraits>::read()
   {
      // Read the Physical parameters
      DynamoFloat E, q, Ra, Ro;
      this->readPhysical(E, q, Ra, Ro);

      // Read the truncation information
      this->readTruncation();

      // Check file compatibility with data truncation
      this->checkTruncation();

      // Set Read arguments
      this->setReadArguments();

      // Read the run information
      this->readRun();

      // Read codensity coefficients
      if(TSimTraits<TSimType>::NeedCodensity)
      {
         this->readCodensity(this->mpCodC->rOc().rPerturbation().data());
      }

      // Read magnetic coefficients
      if(TSimTraits<TSimType>::NeedMagnetic)
      {
         this->readMagnetic(this->mpMagB->rOc().rPerturbation().rTor().data(), this->mpMagB->rOc().rPerturbation().rPol().data());
      }
      
      // Read velocity coefficients
      if(TSimTraits<TSimType>::NeedVelocity)
      {
         this->readVelocity(this->mpVelV->rOc().rPerturbation().rTor().data(), this->mpVelV->rOc().rPerturbation().rPol().data());
      }
   }

}

#endif // STATEFILEREADER_HPP
