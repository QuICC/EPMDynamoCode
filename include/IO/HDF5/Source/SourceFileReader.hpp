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
#include "IO/HDF5/Source/SourceCodMagVelReader.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of HDF5 source field reader
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class SourceFileReader: public SourceCodMagVelReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasVelICSource, TSimTraits<TSimType>::HasVelOCSource>
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

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits> SourceFileReader<TSimType, TSimTraits>::SourceFileReader(typename TSimTraits<TSimType>::CodType &codC, typename TSimTraits<TSimType>::MagType &magB, typename TSimTraits<TSimType>::VelType &velV)
      : SourceCodMagVelReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasVelICSource, TSimTraits<TSimType>::HasVelOCSource>("full", codC.oc().trunc())
   {
      this->mpCodC = &codC;
      // make sure extra values are zero
      this->mpCodC->rOcSrc().initialiseZeros();

      this->mpMagB = &magB;
      this->mpMagB->rOcSrc().initialiseZeros();

      this->mpVelV = &velV;
      this->mpVelV->rOcSrc().initialiseZeros();
   }

   template <typename TSimType, template <typename> class TSimTraits> SourceFileReader<TSimType, TSimTraits>::SourceFileReader(typename TSimTraits<TSimType>::CodType &codC, typename TSimTraits<TSimType>::VelType &velV)
      : SourceCodMagVelReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasVelICSource, TSimTraits<TSimType>::HasVelOCSource>("codVel", codC.oc().trunc())
   {
      this->mpCodC = &codC;
      this->mpCodC->rOcSrc().initialiseZeros();

      this->mpVelV = &velV;
      this->mpVelV->rOcSrc().initialiseZeros();
   }

   template <typename TSimType, template <typename> class TSimTraits> SourceFileReader<TSimType, TSimTraits>::SourceFileReader(typename TSimTraits<TSimType>::MagType &magB, typename TSimTraits<TSimType>::VelType &velV)
      : SourceCodMagVelReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasVelICSource, TSimTraits<TSimType>::HasVelOCSource>("magVel", magB.oc().trunc())
   {
      this->mpMagB = &magB;
      this->mpMagB->rOcSrc().initialiseZeros();

      this->mpVelV = &velV;
      this->mpVelV->rOcSrc().initialiseZeros();
   }

   template <typename TSimType, template <typename> class TSimTraits> SourceFileReader<TSimType, TSimTraits>::SourceFileReader(typename TSimTraits<TSimType>::CodType &codC)
      : SourceCodMagVelReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasVelICSource, TSimTraits<TSimType>::HasVelOCSource>("cod", codC.oc().trunc())
   {
      this->mpCodC = &codC;
      this->mpCodC->rOcSrc().initialiseZeros();
   }

   template <typename TSimType, template <typename> class TSimTraits> SourceFileReader<TSimType, TSimTraits>::SourceFileReader(typename TSimTraits<TSimType>::MagType &magB)
      : SourceCodMagVelReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasVelICSource, TSimTraits<TSimType>::HasVelOCSource>("mag", magB.oc().trunc())
   {
      this->mpMagB = &magB;
      this->mpMagB->rOcSrc().initialiseZeros();
   }

   template <typename TSimType, template <typename> class TSimTraits> SourceFileReader<TSimType, TSimTraits>::SourceFileReader(typename TSimTraits<TSimType>::VelType &velV)
      : SourceCodMagVelReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasVelICSource, TSimTraits<TSimType>::HasVelOCSource>("vel", velV.oc().trunc())
   {
      this->mpVelV = &velV;
      this->mpVelV->rOcSrc().initialiseZeros();
   }

   template <typename TSimType, template <typename> class TSimTraits> void SourceFileReader<TSimType, TSimTraits>::read()
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
