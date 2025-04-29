/** \file ImposedFieldWriter.hpp
 *  \brief Implementation of the HDF5 imposed field writer
 */

#ifndef IMPOSEDFIELDWRITER_HPP
#define IMPOSEDFIELDWRITER_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "IO/HDF5/Imposed/ImposedCodMagVelWriter.hpp"
#include "Timestepping/TimestepParameters.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the HDF5 imposed field writer
    *
    * \tparam TSimTraits Traits of the implementation
    */
   template <typename TSimTraits> class ImposedFieldWriter: public ImposedCodMagVelWriter<TSimTraits, TSimTraits::HasVelICImposed, TSimTraits::HasVelOCImposed>
   {
      public:
         /**
         * @brief Constructor with all three fields
         *
         * @param codC Codensity variable
         * @param magB Magnetic variable
         * @param velV Velocity variable
         */
         ImposedFieldWriter(const typename TSimTraits::CodType &codC, const typename TSimTraits::MagType &magB, const typename TSimTraits::VelType &velV);

         /**
         * @brief Constructor for thermal convection
         *
         * @param codC Codensity variable
         * @param velV Velocity variable
         */
         ImposedFieldWriter(const typename TSimTraits::CodType &codC, const typename TSimTraits::VelType &velV);

         /**
         * @brief Constructor for magneto-convection
         *
         * @param magB Magnetic variable
         * @param velV Velocity variable
         */
         ImposedFieldWriter(const typename TSimTraits::MagType &magB, const typename TSimTraits::VelType &velV);

         /**
         * @brief Constructor with only codensity
         *
         * @param codC Codensity variable
         */
         ImposedFieldWriter(const typename TSimTraits::CodType &codC);

         /**
         * @brief Constructor with only magnetic field
         *
         * @param magB Magnetic variable
         */
         ImposedFieldWriter(const typename TSimTraits::MagType &magB);

         /**
         * @brief Constructor for cases with all three fields
         *
         * @param velV Velocity variable
         */
         ImposedFieldWriter(const typename TSimTraits::VelType &velV);

         /**
         * @brief Destructor
         */
         virtual ~ImposedFieldWriter() {};

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

   template <typename TSimTraits> ImposedFieldWriter<TSimTraits>::ImposedFieldWriter(const typename TSimTraits::CodType &codC, const typename TSimTraits::MagType &magB, const typename TSimTraits::VelType &velV)
      : ImposedCodMagVelWriter<TSimTraits, TSimTraits::HasVelICImposed, TSimTraits::HasVelOCImposed>("full", codC.oc().trunc()), mTrunc(4)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();

      this->mpCodC = &codC;

      this->mpMagB = &magB;

      this->mpVelV = &velV;
   }

   template <typename TSimTraits> ImposedFieldWriter<TSimTraits>::ImposedFieldWriter(const typename TSimTraits::CodType &codC, const typename TSimTraits::VelType &velV)
      : ImposedCodMagVelWriter<TSimTraits, TSimTraits::HasVelICImposed, TSimTraits::HasVelOCImposed>("codVel", codC.oc().trunc()), mTrunc(4)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();

      this->mpCodC = &codC;

      this->mpVelV = &velV;
   }

   template <typename TSimTraits> ImposedFieldWriter<TSimTraits>::ImposedFieldWriter(const typename TSimTraits::MagType &magB, const typename TSimTraits::VelType &velV)
      : ImposedCodMagVelWriter<TSimTraits, TSimTraits::HasVelICImposed, TSimTraits::HasVelOCImposed>("magVel", magB.oc().trunc()), mTrunc(4)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();

      this->mpMagB = &magB;

      this->mpVelV = &velV;
   }

   template <typename TSimTraits> ImposedFieldWriter<TSimTraits>::ImposedFieldWriter(const typename TSimTraits::CodType &codC)
      : ImposedCodMagVelWriter<TSimTraits, TSimTraits::HasVelICImposed, TSimTraits::HasVelOCImposed>("cod", codC.oc().trunc()), mTrunc(4)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();

      this->mpCodC = &codC;
   }

   template <typename TSimTraits> ImposedFieldWriter<TSimTraits>::ImposedFieldWriter(const typename TSimTraits::MagType &magB)
      : ImposedCodMagVelWriter<TSimTraits, TSimTraits::HasVelICImposed, TSimTraits::HasVelOCImposed>("mag", magB.oc().trunc()), mTrunc(4)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();

      this->mpMagB = &magB;
   }

   template <typename TSimTraits> ImposedFieldWriter<TSimTraits>::ImposedFieldWriter(const typename TSimTraits::VelType &velV)
      : ImposedCodMagVelWriter<TSimTraits, TSimTraits::HasVelICImposed, TSimTraits::HasVelOCImposed>("vel", velV.oc().trunc()), mTrunc(4)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();

      this->mpVelV = &velV;
   }

   template <typename TSimTraits> void ImposedFieldWriter<TSimTraits>::write()
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

#endif // IMPOSEDFIELDWRITER_HPP
