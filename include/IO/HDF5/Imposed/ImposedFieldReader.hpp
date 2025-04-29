/** \file ImposedFieldReader.hpp
 *  \brief Implementation of HDF5 imposed field reader
 */

#ifndef IMPOSEDFIELDREADER_HPP
#define IMPOSEDFIELDREADER_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "IO/HDF5/Imposed/ImposedCodMagVelReader.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of HDF5 imposed field reader
    *
    * \tparam TSimTraits Traits of the implementation
    */
   template <typename TSimTraits> class ImposedFieldReader: public ImposedCodMagVelReader<TSimTraits, TSimTraits::HasVelICImposed, TSimTraits::HasVelOCImposed>
   {
      public:
         /**
         * @brief Constructor with all three fields
         *
         * @param codC Codensity variable
         * @param magB Magnetic variable
         * @param velV Velocity variable
         */
         ImposedFieldReader(typename TSimTraits::CodType& codC, typename TSimTraits::MagType& magB, typename TSimTraits::VelType& velV);

         /**
         * @brief Constructor for thermal convection
         *
         * @param codC Codensity variable
         * @param velV Velocity variable
         */
         ImposedFieldReader(typename TSimTraits::CodType& codC, typename TSimTraits::VelType& velV);

         /**
         * @brief Constructor for magneto-convection
         *
         * @param magB Magnetic variable
         * @param velV Velocity variable
         */
         ImposedFieldReader(typename TSimTraits::MagType& magB, typename TSimTraits::VelType& velV);

         /**
         * @brief Constructor with only codensity
         *
         * @param codC Codensity variable
         */
         ImposedFieldReader(typename TSimTraits::CodType& codC);

         /**
         * @brief Constructor with only magnetic field
         *
         * @param magB Magnetic variable
         */
         ImposedFieldReader(typename TSimTraits::MagType& magB);

         /**
         * @brief Constructor for cases with all three fields
         *
         * @param velV Velocity variable
         */
         ImposedFieldReader(typename TSimTraits::VelType& velV);

         /**
         * @brief Destructor
         */
         virtual ~ImposedFieldReader() {};

         /**
          * @brief Read State to file
          */
         virtual void read();
         
      protected:

      private:
   };

   template <typename TSimTraits> ImposedFieldReader<TSimTraits>::ImposedFieldReader(typename TSimTraits::CodType &codC, typename TSimTraits::MagType &magB, typename TSimTraits::VelType &velV)
      : ImposedCodMagVelReader<TSimTraits, TSimTraits::HasVelICImposed, TSimTraits::HasVelOCImposed>("full", codC.oc().trunc())
   {
      this->mpCodC = &codC;

      this->mpMagB = &magB;

      this->mpVelV = &velV;
   }

   template <typename TSimTraits> ImposedFieldReader<TSimTraits>::ImposedFieldReader(typename TSimTraits::CodType &codC, typename TSimTraits::VelType &velV)
      : ImposedCodMagVelReader<TSimTraits, TSimTraits::HasVelICImposed, TSimTraits::HasVelOCImposed>("codVel", codC.oc().trunc())
   {
      this->mpCodC = &codC;

      this->mpVelV = &velV;
   }

   template <typename TSimTraits> ImposedFieldReader<TSimTraits>::ImposedFieldReader(typename TSimTraits::MagType &magB, typename TSimTraits::VelType &velV)
      : ImposedCodMagVelReader<TSimTraits, TSimTraits::HasVelICImposed, TSimTraits::HasVelOCImposed>("magVel", magB.oc().trunc())
   {
      this->mpMagB = &magB;

      this->mpVelV = &velV;
   }

   template <typename TSimTraits> ImposedFieldReader<TSimTraits>::ImposedFieldReader(typename TSimTraits::CodType &codC)
      : ImposedCodMagVelReader<TSimTraits, TSimTraits::HasVelICImposed, TSimTraits::HasVelOCImposed>("cod", codC.oc().trunc())
   {
      this->mpCodC = &codC;
   }

   template <typename TSimTraits> ImposedFieldReader<TSimTraits>::ImposedFieldReader(typename TSimTraits::MagType &magB)
      : ImposedCodMagVelReader<TSimTraits, TSimTraits::HasVelICImposed, TSimTraits::HasVelOCImposed>("mag", magB.oc().trunc())
   {
      this->mpMagB = &magB;
   }

   template <typename TSimTraits> ImposedFieldReader<TSimTraits>::ImposedFieldReader(typename TSimTraits::VelType &velV)
      : ImposedCodMagVelReader<TSimTraits, TSimTraits::HasVelICImposed, TSimTraits::HasVelOCImposed>("vel", velV.oc().trunc())
   {
      this->mpVelV = &velV;
   }

   template <typename TSimTraits> void ImposedFieldReader<TSimTraits>::read()
   {
      // Read the truncation information
      this->readTruncation();

      // Check file compatibility with data truncation
      this->checkTruncation();

      // Set Read arguments
      this->setReadArguments();

      // Read codensity coefficients
      this->readCod();

      // Read magnetic coefficients
      this->readMag();
      
      // Read velocity coefficients
      this->readVel();
   }

}

#endif // IMPOSEDFIELDREADER_HPP
