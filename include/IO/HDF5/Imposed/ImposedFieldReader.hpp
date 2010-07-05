/** \file ImposedFieldReader.hpp
 *  \brief Implementation of HDF5 imposed field reader
 */

#ifndef IMPOSEDFIELDREADER_HPP
#define IMPOSEDFIELDREADER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "IO/HDF5/Imposed/ImposedFieldReaderBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of HDF5 imposed field reader
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class ImposedFieldReader: public ImposedFieldReaderBase
   {
      public:
         /**
         * @brief Constructor with all three fields
         *
         * @param codC Codensity variable
         * @param magB Magnetic variable
         * @param velV Velocity variable
         */
         ImposedFieldReader(typename TSimTraits<TSimType>::CodType& codC, typename TSimTraits<TSimType>::MagType& magB, typename TSimTraits<TSimType>::VelType& velV);

         /**
         * @brief Constructor for thermal convection
         *
         * @param codC Codensity variable
         * @param velV Velocity variable
         */
         ImposedFieldReader(typename TSimTraits<TSimType>::CodType& codC, typename TSimTraits<TSimType>::VelType& velV);

         /**
         * @brief Constructor for magneto-convection
         *
         * @param magB Magnetic variable
         * @param velV Velocity variable
         */
         ImposedFieldReader(typename TSimTraits<TSimType>::MagType& magB, typename TSimTraits<TSimType>::VelType& velV);

         /**
         * @brief Constructor with only codensity
         *
         * @param codC Codensity variable
         */
         ImposedFieldReader(typename TSimTraits<TSimType>::CodType& codC);

         /**
         * @brief Constructor with only magnetic field
         *
         * @param magB Magnetic variable
         */
         ImposedFieldReader(typename TSimTraits<TSimType>::MagType& magB);

         /**
         * @brief Constructor for cases with all three fields
         *
         * @param velV Velocity variable
         */
         ImposedFieldReader(typename TSimTraits<TSimType>::VelType& velV);

         /**
         * @brief Destructor
         */
         virtual ~ImposedFieldReader() {};

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

   template <typename TSimType, template <typename> class TSimTraits> ImposedFieldReader<TSimType, TSimTraits>::ImposedFieldReader(typename TSimTraits<TSimType>::CodType &codC, typename TSimTraits<TSimType>::MagType &magB, typename TSimTraits<TSimType>::VelType &velV)
      : ImposedFieldReaderBase("full", codC.trunc()), mpCodC(&codC), mpMagB(&magB), mpVelV(&velV)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> ImposedFieldReader<TSimType, TSimTraits>::ImposedFieldReader(typename TSimTraits<TSimType>::CodType &codC, typename TSimTraits<TSimType>::VelType &velV)
      : ImposedFieldReaderBase("codVel", codC.trunc()), mpCodC(&codC), mpMagB(NULL), mpVelV(&velV)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> ImposedFieldReader<TSimType, TSimTraits>::ImposedFieldReader(typename TSimTraits<TSimType>::MagType &magB, typename TSimTraits<TSimType>::VelType &velV)
      : ImposedFieldReaderBase("magVel", magB.trunc()), mpCodC(NULL), mpMagB(&magB), mpVelV(&velV)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> ImposedFieldReader<TSimType, TSimTraits>::ImposedFieldReader(typename TSimTraits<TSimType>::CodType &codC)
      : ImposedFieldReaderBase("cod", codC.trunc()), mpCodC(&codC), mpMagB(NULL), mpVelV(NULL)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> ImposedFieldReader<TSimType, TSimTraits>::ImposedFieldReader(typename TSimTraits<TSimType>::MagType &magB)
      : ImposedFieldReaderBase("mag", magB.trunc()), mpCodC(NULL), mpMagB(&magB), mpVelV(NULL)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> ImposedFieldReader<TSimType, TSimTraits>::ImposedFieldReader(typename TSimTraits<TSimType>::VelType &velV)
      : ImposedFieldReaderBase("vel", velV.trunc()), mpCodC(NULL), mpMagB(NULL), mpVelV(&velV)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void ImposedFieldReader<TSimType, TSimTraits>::read()
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
         this->readCodensity(this->mpCodC->rOc().rImposed().data());
      }

      // Read magnetic coefficients
      if(this->mpMagB != NULL)
      {
         this->readMagnetic(this->mpMagB->rOc().rImposed().rTor().data(), this->mpMagB->rOc().rImposed().rPol().data());
      }
      
      // Read velocity coefficients
      if(this->mpVelV != NULL)
      {
         this->readVelocity(this->mpVelV->rOc().rImposed().rTor().data(), this->mpVelV->rOc().rImposed().rPol().data());
      }
   }

}

#endif // IMPOSEDFIELDREADER_HPP
