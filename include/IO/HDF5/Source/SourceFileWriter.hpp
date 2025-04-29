/** \file SourceFileWriter.hpp
 *  \brief Implementation of the HDF5 source field writer
 */

#ifndef SOURCEFILEWRITER_HPP
#define SOURCEFILEWRITER_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "IO/HDF5/Source/SourceCodMagVelWriter.hpp"
#include "Timestepping/TimestepParameters.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the HDF5 source field writer
    *
    * \tparam TSimTraits Traits of the implementation
    */
   template <typename TSimTraits> class SourceFileWriter: public SourceCodMagVelWriter<TSimTraits, TSimTraits::HasVelICSource, TSimTraits::HasVelOCSource>
   {
      public:
         /**
         * @brief Constructor with all three fields
         *
         * @param codC Codensity variable
         * @param magB Magnetic variable
         * @param velV Velocity variable
         */
         SourceFileWriter(const typename TSimTraits::CodType &codC, const typename TSimTraits::MagType &magB, const typename TSimTraits::VelType &velV);

         /**
         * @brief Constructor for thermal convection
         *
         * @param codC Codensity variable
         * @param velV Velocity variable
         */
         SourceFileWriter(const typename TSimTraits::CodType &codC, const typename TSimTraits::VelType &velV);

         /**
         * @brief Constructor for magneto-convection
         *
         * @param magB Magnetic variable
         * @param velV Velocity variable
         */
         SourceFileWriter(const typename TSimTraits::MagType &magB, const typename TSimTraits::VelType &velV);

         /**
         * @brief Constructor with only codensity
         *
         * @param codC Codensity variable
         */
         SourceFileWriter(const typename TSimTraits::CodType &codC);

         /**
         * @brief Constructor with only magnetic field
         *
         * @param magB Magnetic variable
         */
         SourceFileWriter(const typename TSimTraits::MagType &magB);

         /**
         * @brief Constructor for cases with all three fields
         *
         * @param velV Velocity variable
         */
         SourceFileWriter(const typename TSimTraits::VelType &velV);

         /**
         * @brief Destructor
         */
         virtual ~SourceFileWriter() {};

         /**
          * @brief Write State to file
          */
         virtual void write();
         
      protected:

      private:
         /**
          * @brief Store truncation information in a array for faster access
          */
         ArrayI   mTrunc;
   };

   template <typename TSimTraits> SourceFileWriter<TSimTraits>::SourceFileWriter(const typename TSimTraits::CodType &codC, const typename TSimTraits::MagType &magB, const typename TSimTraits::VelType &velV)
      : SourceCodMagVelWriter<TSimTraits, TSimTraits::HasVelICSource, TSimTraits::HasVelOCSource>("full", codC.oc().trunc()), mTrunc(4)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();

      this->mpCodC = &codC;

      this->mpMagB = &magB;

      this->mpVelV = &velV;
   }

   template <typename TSimTraits> SourceFileWriter<TSimTraits>::SourceFileWriter(const typename TSimTraits::CodType &codC, const typename TSimTraits::VelType &velV)
      : SourceCodMagVelWriter<TSimTraits, TSimTraits::HasVelICSource, TSimTraits::HasVelOCSource>("codVel", codC.oc().trunc()), mTrunc(4)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();

      this->mpCodC = &codC;

      this->mpVelV = &velV;
   }

   template <typename TSimTraits> SourceFileWriter<TSimTraits>::SourceFileWriter(const typename TSimTraits::MagType &magB, const typename TSimTraits::VelType &velV)
      : SourceCodMagVelWriter<TSimTraits, TSimTraits::HasVelICSource, TSimTraits::HasVelOCSource>("magVel", magB.oc().trunc()), mTrunc(4)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();

      this->mpMagB = &magB;

      this->mpVelV = &velV;
   }

   template <typename TSimTraits> SourceFileWriter<TSimTraits>::SourceFileWriter(const typename TSimTraits::CodType &codC)
      : SourceCodMagVelWriter<TSimTraits, TSimTraits::HasVelICSource, TSimTraits::HasVelOCSource>("cod", codC.oc().trunc()), mTrunc(4)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();

      this->mpCodC = &codC;
   }

   template <typename TSimTraits> SourceFileWriter<TSimTraits>::SourceFileWriter(const typename TSimTraits::MagType &magB)
      : SourceCodMagVelWriter<TSimTraits, TSimTraits::HasVelICSource, TSimTraits::HasVelOCSource>("mag", magB.oc().trunc()), mTrunc(4)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();

      this->mpMagB = &magB;
   }

   template <typename TSimTraits> SourceFileWriter<TSimTraits>::SourceFileWriter(const typename TSimTraits::VelType &velV)
      : SourceCodMagVelWriter<TSimTraits, TSimTraits::HasVelICSource, TSimTraits::HasVelOCSource>("vel", velV.oc().trunc()), mTrunc(4)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();

      this->mpVelV = &velV;
   }

   template <typename TSimTraits> void SourceFileWriter<TSimTraits>::write()
   {
      // Create file
      this->preWrite();

      // Create the header and version information
      this->createFileInfo();

      // Write the truncation information
      this->writeTruncation(this->mTrunc(0), this->mTrunc(1), this->mTrunc(2), this->mTrunc(3));

      // Write the codensity coefficients
      this->writeCod();

      // Write the magnetic coefficients
      this->writeMag();

      // Write the velocity coefficients
      this->writeVel();

      // Close file
      this->postWrite();
   }

}

#endif // SOURCEFILEWRITER_HPP
