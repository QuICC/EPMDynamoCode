/** \file SourceFileReader.hpp
 *  \brief Implementation of HDF5 source field reader
 */

#ifndef SOURCEFILEREADER_HPP
#define SOURCEFILEREADER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "IO/HDF5/Source/SourceFileReaderBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of HDF5 source field reader
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class SourceFileReader: public SourceFileReaderBase
   {
      public:
         /**
         * @brief Constructor with all three fields
         *
         * @param codC Codensity variable
         * @param magB Magnetic variable
         * @param velV Velocity variable
         */
         SourceFileReader(typename TSimTraits<TSimType>::CodType& codC, typename TSimTraits<TSimType>::MagType& magB, typename TSimTraits<TSimType>::VelType& velV);

         /**
         * @brief Constructor for thermal convection
         *
         * @param codC Codensity variable
         * @param velV Velocity variable
         */
         SourceFileReader(typename TSimTraits<TSimType>::CodType& codC, typename TSimTraits<TSimType>::VelType& velV);

         /**
         * @brief Constructor for magneto-convection
         *
         * @param magB Magnetic variable
         * @param velV Velocity variable
         */
         SourceFileReader(typename TSimTraits<TSimType>::MagType& magB, typename TSimTraits<TSimType>::VelType& velV);

         /**
         * @brief Constructor with only codensity
         *
         * @param codC Codensity variable
         */
         SourceFileReader(typename TSimTraits<TSimType>::CodType& codC);

         /**
         * @brief Constructor with only magnetic field
         *
         * @param magB Magnetic variable
         */
         SourceFileReader(typename TSimTraits<TSimType>::MagType& magB);

         /**
         * @brief Constructor for cases with all three fields
         *
         * @param velV Velocity variable
         */
         SourceFileReader(typename TSimTraits<TSimType>::VelType& velV);

         /**
         * @brief Destructor
         */
         virtual ~SourceFileReader() {};

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

   template <typename TSimType, template <typename> class TSimTraits> SourceFileReader<TSimType, TSimTraits>::SourceFileReader(typename TSimTraits<TSimType>::CodType &codC, typename TSimTraits<TSimType>::MagType &magB, typename TSimTraits<TSimType>::VelType &velV)
      : SourceFileReaderBase("full", codC.trunc()), mpCodC(&codC), mpMagB(&magB), mpVelV(&velV)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> SourceFileReader<TSimType, TSimTraits>::SourceFileReader(typename TSimTraits<TSimType>::CodType &codC, typename TSimTraits<TSimType>::VelType &velV)
      : SourceFileReaderBase("codVel", codC.trunc()), mpCodC(&codC), mpMagB(NULL), mpVelV(&velV)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> SourceFileReader<TSimType, TSimTraits>::SourceFileReader(typename TSimTraits<TSimType>::MagType &magB, typename TSimTraits<TSimType>::VelType &velV)
      : SourceFileReaderBase("magVel", magB.trunc()), mpCodC(NULL), mpMagB(&magB), mpVelV(&velV)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> SourceFileReader<TSimType, TSimTraits>::SourceFileReader(typename TSimTraits<TSimType>::CodType &codC)
      : SourceFileReaderBase("cod", codC.trunc()), mpCodC(&codC), mpMagB(NULL), mpVelV(NULL)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> SourceFileReader<TSimType, TSimTraits>::SourceFileReader(typename TSimTraits<TSimType>::MagType &magB)
      : SourceFileReaderBase("mag", magB.trunc()), mpCodC(NULL), mpMagB(&magB), mpVelV(NULL)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> SourceFileReader<TSimType, TSimTraits>::SourceFileReader(typename TSimTraits<TSimType>::VelType &velV)
      : SourceFileReaderBase("vel", velV.trunc()), mpCodC(NULL), mpMagB(NULL), mpVelV(&velV)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void SourceFileReader<TSimType, TSimTraits>::read()
   {
      // Read the truncation information
      this->readTruncation();

      // Check file compatibility with data truncation
      this->checkTruncation();

      // Set Read arguments
      this->setReadArguments();

      // Read codensity coefficients
      if(this->mpCodC != NULL)
      {
         this->readCodensity(this->mpCodC->rOc().rPerturbation().data());
      }

      // Read magnetic coefficients
      if(this->mpMagB != NULL)
      {
         this->readMagnetic(this->mpMagB->rOc().rPerturbation().rTor().data(), this->mpMagB->rOc().rPerturbation().rPol().data());
      }
      
      // Read velocity coefficients
      if(this->mpVelV != NULL)
      {
         this->readVelocity(this->mpVelV->rOc().rPerturbation().rTor().data(), this->mpVelV->rOc().rPerturbation().rPol().data());
      }
   }

}

#endif // SOURCEFILEREADER_HPP
