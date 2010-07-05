/** \file ImposedFieldWriter.hpp
 *  \brief Implementation of the HDF5 imposed field writer
 */

#ifndef IMPOSEDFIELDWRITER_HPP
#define IMPOSEDFIELDWRITER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "IO/HDF5/Imposed/ImposedFieldWriterBase.hpp"
#include "Timestepping/TimestepParameters.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the HDF5 imposed field writer
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class ImposedFieldWriter: public ImposedFieldWriterBase
   {
      public:
         /**
         * @brief Constructor with all three fields
         *
         * @param codC Codensity variable
         * @param magB Magnetic variable
         * @param velV Velocity variable
         */
         ImposedFieldWriter(const typename TSimTraits<TSimType>::CodType &codC, const typename TSimTraits<TSimType>::MagType &magB, const typename TSimTraits<TSimType>::VelType &velV);

         /**
         * @brief Constructor for thermal convection
         *
         * @param codC Codensity variable
         * @param velV Velocity variable
         */
         ImposedFieldWriter(const typename TSimTraits<TSimType>::CodType &codC, const typename TSimTraits<TSimType>::VelType &velV);

         /**
         * @brief Constructor for magneto-convection
         *
         * @param magB Magnetic variable
         * @param velV Velocity variable
         */
         ImposedFieldWriter(const typename TSimTraits<TSimType>::MagType &magB, const typename TSimTraits<TSimType>::VelType &velV);

         /**
         * @brief Constructor with only codensity
         *
         * @param codC Codensity variable
         */
         ImposedFieldWriter(const typename TSimTraits<TSimType>::CodType &codC);

         /**
         * @brief Constructor with only magnetic field
         *
         * @param magB Magnetic variable
         */
         ImposedFieldWriter(const typename TSimTraits<TSimType>::MagType &magB);

         /**
         * @brief Constructor for cases with all three fields
         *
         * @param velV Velocity variable
         */
         ImposedFieldWriter(const typename TSimTraits<TSimType>::VelType &velV);

         /**
         * @brief Destructor
         */
         virtual ~ImposedFieldWriter() {};

         /**
          * @brief Write State to file
          */
         virtual void write();
         
      protected:
         /**
          * @brief Pointer to the condensity variable
          */
         const typename TSimTraits<TSimType>::CodType*  mpCodC;

         /**
          * @brief Pointer to the magnetic variable
          */
         const typename TSimTraits<TSimType>::MagType*  mpMagB;

         /**
          * @brief Pointer to the velocity variable
          */
         const typename TSimTraits<TSimType>::VelType*  mpVelV;

      private:
         /**
          * @brief Store truncation information in a array for faster access
          */
         ArrayI   mTrunc;
   };

   template <typename TSimType, template <typename> class TSimTraits> ImposedFieldWriter<TSimType, TSimTraits>::ImposedFieldWriter(const typename TSimTraits<TSimType>::CodType &codC, const typename TSimTraits<TSimType>::MagType &magB, const typename TSimTraits<TSimType>::VelType &velV)
      : ImposedFieldWriterBase("full", codC.trunc()), mpCodC(&codC), mpMagB(&magB), mpVelV(&velV), mTrunc(4)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();
   }

   template <typename TSimType, template <typename> class TSimTraits> ImposedFieldWriter<TSimType, TSimTraits>::ImposedFieldWriter(const typename TSimTraits<TSimType>::CodType &codC, const typename TSimTraits<TSimType>::VelType &velV)
      : ImposedFieldWriterBase("codVel", codC.trunc()), mpCodC(&codC), mpMagB(NULL), mpVelV(&velV), mTrunc(4)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();
   }

   template <typename TSimType, template <typename> class TSimTraits> ImposedFieldWriter<TSimType, TSimTraits>::ImposedFieldWriter(const typename TSimTraits<TSimType>::MagType &magB, const typename TSimTraits<TSimType>::VelType &velV)
      : ImposedFieldWriterBase("magVel", magB.trunc()), mpCodC(NULL), mpMagB(&magB), mpVelV(&velV), mTrunc(4)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();
   }

   template <typename TSimType, template <typename> class TSimTraits> ImposedFieldWriter<TSimType, TSimTraits>::ImposedFieldWriter(const typename TSimTraits<TSimType>::CodType &codC)
      : ImposedFieldWriterBase("cod", codC.trunc()), mpCodC(&codC), mpMagB(NULL), mpVelV(NULL), mTrunc(4)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();
   }

   template <typename TSimType, template <typename> class TSimTraits> ImposedFieldWriter<TSimType, TSimTraits>::ImposedFieldWriter(const typename TSimTraits<TSimType>::MagType &magB)
      : ImposedFieldWriterBase("mag", magB.trunc()), mpCodC(NULL), mpMagB(&magB), mpVelV(NULL), mTrunc(4)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();
   }

   template <typename TSimType, template <typename> class TSimTraits> ImposedFieldWriter<TSimType, TSimTraits>::ImposedFieldWriter(const typename TSimTraits<TSimType>::VelType &velV)
      : ImposedFieldWriterBase("vel", velV.trunc()), mpCodC(NULL), mpMagB(NULL), mpVelV(&velV), mTrunc(4)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();
   }

   template <typename TSimType, template <typename> class TSimTraits> void ImposedFieldWriter<TSimType, TSimTraits>::write()
   {
      // Create file
      this->preWrite();

      // Create the header and version information
      this->createFileInfo();

      // Write the truncation information
      this->writeTruncation(this->mTrunc(0), this->mTrunc(1), this->mTrunc(2), this->mTrunc(3));

      // Write the codensity coefficients
      if(this->mpCodC != NULL)
      {
         this->writeCodensity(this->mpCodC->oc().imposed().data());
      }

      // Write the magnetic coefficients
      if(this->mpMagB != NULL)
      {
         this->writeMagnetic(this->mpMagB->oc().imposed().tor().data(), this->mpMagB->oc().imposed().pol().data());
      }

      // Write the velocity coefficients
      if(this->mpVelV != NULL)
      {
         this->writeVelocity(this->mpVelV->oc().imposed().tor().data(), this->mpVelV->oc().imposed().pol().data());
      }

      // Close file
      this->postWrite();
   }

}

#endif // IMPOSEDFIELDWRITER_HPP
