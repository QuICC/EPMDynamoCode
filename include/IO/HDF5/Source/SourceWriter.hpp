/** \file SourceWriter.hpp
 *  \brief Implementation of the HDF5 source field writer
 */

#ifndef SOURCEWRITER_HPP
#define SOURCEWRITER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "IO/HDF5/Source/SourceWriterBase.hpp"
#include "Timestepping/TimestepParameters.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the HDF5 source field writer
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class SourceWriter: public SourceWriterBase
   {
      public:
         /**
         * @brief Constructor with all three fields
         *
         * @param codC Codensity variable
         * @param magB Magnetic variable
         * @param velV Velocity variable
         */
         SourceWriter(const typename TSimTraits<TSimType>::CodType &codC, const typename TSimTraits<TSimType>::MagType &magB, const typename TSimTraits<TSimType>::VelType &velV);

         /**
         * @brief Constructor for thermal convection
         *
         * @param codC Codensity variable
         * @param velV Velocity variable
         */
         SourceWriter(const typename TSimTraits<TSimType>::CodType &codC, const typename TSimTraits<TSimType>::VelType &velV);

         /**
         * @brief Constructor for magneto-convection
         *
         * @param magB Magnetic variable
         * @param velV Velocity variable
         */
         SourceWriter(const typename TSimTraits<TSimType>::MagType &magB, const typename TSimTraits<TSimType>::VelType &velV);

         /**
         * @brief Constructor with only codensity
         *
         * @param codC Codensity variable
         */
         SourceWriter(const typename TSimTraits<TSimType>::CodType &codC);

         /**
         * @brief Constructor with only magnetic field
         *
         * @param magB Magnetic variable
         */
         SourceWriter(const typename TSimTraits<TSimType>::MagType &magB);

         /**
         * @brief Constructor for cases with all three fields
         *
         * @param velV Velocity variable
         */
         SourceWriter(const typename TSimTraits<TSimType>::VelType &velV);

         /**
         * @brief Destructor
         */
         virtual ~SourceWriter() {};

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

   template <typename TSimType, template <typename> class TSimTraits> SourceWriter<TSimType, TSimTraits>::SourceWriter(const typename TSimTraits<TSimType>::CodType &codC, const typename TSimTraits<TSimType>::MagType &magB, const typename TSimTraits<TSimType>::VelType &velV)
      : SourceWriterBase("full", codC.trunc()), mpCodC(&codC), mpMagB(&magB), mpVelV(&velV), mTrunc(4)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();
   }

   template <typename TSimType, template <typename> class TSimTraits> SourceWriter<TSimType, TSimTraits>::SourceWriter(const typename TSimTraits<TSimType>::CodType &codC, const typename TSimTraits<TSimType>::VelType &velV)
      : SourceWriterBase("codVel", codC.trunc()), mpCodC(&codC), mpMagB(NULL), mpVelV(&velV), mTrunc(4)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();
   }

   template <typename TSimType, template <typename> class TSimTraits> SourceWriter<TSimType, TSimTraits>::SourceWriter(const typename TSimTraits<TSimType>::MagType &magB, const typename TSimTraits<TSimType>::VelType &velV)
      : SourceWriterBase("magVel", magB.trunc()), mpCodC(NULL), mpMagB(&magB), mpVelV(&velV), mTrunc(4)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();
   }

   template <typename TSimType, template <typename> class TSimTraits> SourceWriter<TSimType, TSimTraits>::SourceWriter(const typename TSimTraits<TSimType>::CodType &codC)
      : SourceWriterBase("cod", codC.trunc()), mpCodC(&codC), mpMagB(NULL), mpVelV(NULL), mTrunc(4)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();
   }

   template <typename TSimType, template <typename> class TSimTraits> SourceWriter<TSimType, TSimTraits>::SourceWriter(const typename TSimTraits<TSimType>::MagType &magB)
      : SourceWriterBase("mag", magB.trunc()), mpCodC(NULL), mpMagB(&magB), mpVelV(NULL), mTrunc(4)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();
   }

   template <typename TSimType, template <typename> class TSimTraits> SourceWriter<TSimType, TSimTraits>::SourceWriter(const typename TSimTraits<TSimType>::VelType &velV)
      : SourceWriterBase("vel", velV.trunc()), mpCodC(NULL), mpMagB(NULL), mpVelV(&velV), mTrunc(4)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();
   }

   template <typename TSimType, template <typename> class TSimTraits> void SourceWriter<TSimType, TSimTraits>::write()
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
         this->writeCodensity(this->mpCodC->oc().perturbation().data());
      }

      // Write the magnetic coefficients
      if(this->mpMagB != NULL)
      {
         this->writeMagnetic(this->mpMagB->oc().perturbation().tor().data(), this->mpMagB->oc().perturbation().pol().data());
      }

      // Write the velocity coefficients
      if(this->mpVelV != NULL)
      {
         this->writeVelocity(this->mpVelV->oc().perturbation().tor().data(), this->mpVelV->oc().perturbation().pol().data());
      }

      // Close file
      this->postWrite();
   }

}

#endif // SOURCEWRITER_HPP
