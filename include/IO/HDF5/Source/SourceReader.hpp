/** \file SourceReader.hpp
 *  \brief Implementation of HDF5 source field reader
 */

#ifndef SOURCEREADER_HPP
#define SOURCEREADER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "IO/HDF5/SourceReaderBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of HDF5 source field reader
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class SourceReader: public SourceReaderBase
   {
      public:
         /**
         * @brief Constructor with all three fields
         *
         * @param codC Codensity variable
         * @param magB Magnetic variable
         * @param velV Velocity variable
         */
         SourceReader(typename TSimTraits<TSimType>::CodType& codC, typename TSimTraits<TSimType>::MagType& magB, typename TSimTraits<TSimType>::VelType& velV);

         /**
         * @brief Constructor for thermal convection
         *
         * @param codC Codensity variable
         * @param velV Velocity variable
         */
         SourceReader(typename TSimTraits<TSimType>::CodType& codC, typename TSimTraits<TSimType>::VelType& velV);

         /**
         * @brief Constructor for magneto-convection
         *
         * @param magB Magnetic variable
         * @param velV Velocity variable
         */
         SourceReader(typename TSimTraits<TSimType>::MagType& magB, typename TSimTraits<TSimType>::VelType& velV);

         /**
         * @brief Constructor with only codensity
         *
         * @param codC Codensity variable
         */
         SourceReader(typename TSimTraits<TSimType>::CodType& codC);

         /**
         * @brief Constructor with only magnetic field
         *
         * @param magB Magnetic variable
         */
         SourceReader(typename TSimTraits<TSimType>::MagType& magB);

         /**
         * @brief Constructor for cases with all three fields
         *
         * @param velV Velocity variable
         */
         SourceReader(typename TSimTraits<TSimType>::VelType& velV);

         /**
         * @brief Destructor
         */
         virtual ~SourceReader() {};

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

   template <typename TSimType, template <typename> class TSimTraits> SourceReader<TSimType, TSimTraits>::SourceReader(typename TSimTraits<TSimType>::CodType &codC, typename TSimTraits<TSimType>::MagType &magB, typename TSimTraits<TSimType>::VelType &velV, std::string name)
      : SourceReaderBase("full", codC.trunc()), mpCodC(codC), mpMagB(magB), mpVelV(velV)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> SourceReader<TSimType, TSimTraits>::SourceReader(typename TSimTraits<TSimType>::CodType &codC, typename TSimTraits<TSimType>::VelType &velV, std::string name)
      : SourceReaderBase("codVel", codC.trunc()), mpCodC(codC), mpMagB(NULL), mpVelV(velV)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> SourceReader<TSimType, TSimTraits>::SourceReader(typename TSimTraits<TSimType>::MagType &magB, typename TSimTraits<TSimType>::VelType &velV, std::string name)
      : SourceReaderBase("magVel", codC.trunc()), mpCodC(NULL), mpMagB(magB), mpVelV(velV)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> SourceReader<TSimType, TSimTraits>::SourceReader(typename TSimTraits<TSimType>::CodType &codC, std::string name)
      : SourceReaderBase("cod", codC.trunc()), mpCodC(codC), mpMagB(NULL), mpVelV(NULL)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> SourceReader<TSimType, TSimTraits>::SourceReader(typename TSimTraits<TSimType>::MagType &magB, std::string name)
      : SourceReaderBase("mag", magB.trunc()), mpCodC(NULL), mpMagB(magB), mpVelV(NULL)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> SourceReader<TSimType, TSimTraits>::SourceReader(typename TSimTraits<TSimType>::VelType &codC, std::string name)
      : SourceReaderBase("vel", velV.trunc()), mpCodC(NULL), mpMagB(NULL), mpVelV(velV)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void SourceReader<TSimType, TSimTraits>::read()
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

#endif // SOURCEREADER_HPP
