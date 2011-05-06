/** \file SourceFileReader.hpp
 *  \brief Implementation of HDF5 source field reader
 */

#ifndef SOURCEFILEREADER_HPP
#define SOURCEFILEREADER_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "IO/HDF5/Source/SourceCodMagVelReader.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of HDF5 source field reader
    *
    * \tparam TSimTraits Traits of the implementation
    */
   template <typename TSimTraits> class SourceFileReader: public SourceCodMagVelReader<TSimTraits, TSimTraits::HasVelICSource, TSimTraits::HasVelOCSource>
   {
      public:
         /**
         * @brief Constructor with all three fields
         *
         * @param codC Codensity variable
         * @param magB Magnetic variable
         * @param velV Velocity variable
         */
         SourceFileReader(typename TSimTraits::CodType& codC, typename TSimTraits::MagType& magB, typename TSimTraits::VelType& velV);

         /**
         * @brief Constructor for thermal convection
         *
         * @param codC Codensity variable
         * @param velV Velocity variable
         */
         SourceFileReader(typename TSimTraits::CodType& codC, typename TSimTraits::VelType& velV);

         /**
         * @brief Constructor for magneto-convection
         *
         * @param magB Magnetic variable
         * @param velV Velocity variable
         */
         SourceFileReader(typename TSimTraits::MagType& magB, typename TSimTraits::VelType& velV);

         /**
         * @brief Constructor with only codensity
         *
         * @param codC Codensity variable
         */
         SourceFileReader(typename TSimTraits::CodType& codC);

         /**
         * @brief Constructor with only magnetic field
         *
         * @param magB Magnetic variable
         */
         SourceFileReader(typename TSimTraits::MagType& magB);

         /**
         * @brief Constructor for cases with all three fields
         *
         * @param velV Velocity variable
         */
         SourceFileReader(typename TSimTraits::VelType& velV);

         /**
         * @brief Destructor
         */
         virtual ~SourceFileReader() {};

         /**
          * @brief Read State to file
          */
         virtual void read();
         
      protected:

      private:
   };

   template <typename TSimTraits> SourceFileReader<TSimTraits>::SourceFileReader(typename TSimTraits::CodType &codC, typename TSimTraits::MagType &magB, typename TSimTraits::VelType &velV)
      : SourceCodMagVelReader<TSimTraits, TSimTraits::HasVelICSource, TSimTraits::HasVelOCSource>("full", codC.oc().trunc())
   {
      this->mpCodC = &codC;
      // make sure extra values are zero
      this->mpCodC->rOcSrc().initialiseZeros();

      this->mpMagB = &magB;
      this->mpMagB->rOcSrc().initialiseZeros();

      this->mpVelV = &velV;
      this->mpVelV->rOcSrc().initialiseZeros();
   }

   template <typename TSimTraits> SourceFileReader<TSimTraits>::SourceFileReader(typename TSimTraits::CodType &codC, typename TSimTraits::VelType &velV)
      : SourceCodMagVelReader<TSimTraits, TSimTraits::HasVelICSource, TSimTraits::HasVelOCSource>("codVel", codC.oc().trunc())
   {
      this->mpCodC = &codC;
      this->mpCodC->rOcSrc().initialiseZeros();

      this->mpVelV = &velV;
      this->mpVelV->rOcSrc().initialiseZeros();
   }

   template <typename TSimTraits> SourceFileReader<TSimTraits>::SourceFileReader(typename TSimTraits::MagType &magB, typename TSimTraits::VelType &velV)
      : SourceCodMagVelReader<TSimTraits, TSimTraits::HasVelICSource, TSimTraits::HasVelOCSource>("magVel", magB.oc().trunc())
   {
      this->mpMagB = &magB;
      this->mpMagB->rOcSrc().initialiseZeros();

      this->mpVelV = &velV;
      this->mpVelV->rOcSrc().initialiseZeros();
   }

   template <typename TSimTraits> SourceFileReader<TSimTraits>::SourceFileReader(typename TSimTraits::CodType &codC)
      : SourceCodMagVelReader<TSimTraits, TSimTraits::HasVelICSource, TSimTraits::HasVelOCSource>("cod", codC.oc().trunc())
   {
      this->mpCodC = &codC;
      this->mpCodC->rOcSrc().initialiseZeros();
   }

   template <typename TSimTraits> SourceFileReader<TSimTraits>::SourceFileReader(typename TSimTraits::MagType &magB)
      : SourceCodMagVelReader<TSimTraits, TSimTraits::HasVelICSource, TSimTraits::HasVelOCSource>("mag", magB.oc().trunc())
   {
      this->mpMagB = &magB;
      this->mpMagB->rOcSrc().initialiseZeros();
   }

   template <typename TSimTraits> SourceFileReader<TSimTraits>::SourceFileReader(typename TSimTraits::VelType &velV)
      : SourceCodMagVelReader<TSimTraits, TSimTraits::HasVelICSource, TSimTraits::HasVelOCSource>("vel", velV.oc().trunc())
   {
      this->mpVelV = &velV;
      this->mpVelV->rOcSrc().initialiseZeros();
   }

   template <typename TSimTraits> void SourceFileReader<TSimTraits>::read()
   {
      // Read the truncation information
      this->readTruncation();

      // Check file compatibility with data truncation
      this->checkTruncation();

      // Set Read arguments
      this->setReadArguments();

      // Read codensity coefficients (if applicable)
      this->readCod();

      // Read magnetic coefficients (if applicable)
      this->readMag();

      // Read velocity coefficients (if applicable)
      this->readVel();
   }

}

#endif // SOURCEFILEREADER_HPP
