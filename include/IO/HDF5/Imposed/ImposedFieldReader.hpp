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
#include "IO/HDF5/Imposed/ImposedCodMagVelReader.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of HDF5 imposed field reader
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class ImposedFieldReader: public ImposedCodMagVelReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasVelICImposed, TSimTraits<TSimType>::HasVelOCImposed>
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

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits> ImposedFieldReader<TSimType, TSimTraits>::ImposedFieldReader(typename TSimTraits<TSimType>::CodType &codC, typename TSimTraits<TSimType>::MagType &magB, typename TSimTraits<TSimType>::VelType &velV)
      : ImposedCodMagVelReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasVelICImposed, TSimTraits<TSimType>::HasVelOCImposed>("full", codC.oc().trunc())
   {
      this->mpCodC = &codC;

      this->mpMagB = &magB;

      this->mpVelV = &velV;
   }

   template <typename TSimType, template <typename> class TSimTraits> ImposedFieldReader<TSimType, TSimTraits>::ImposedFieldReader(typename TSimTraits<TSimType>::CodType &codC, typename TSimTraits<TSimType>::VelType &velV)
      : ImposedCodMagVelReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasVelICImposed, TSimTraits<TSimType>::HasVelOCImposed>("codVel", codC.oc().trunc())
   {
      this->mpCodC = &codC;

      this->mpVelV = &velV;
   }

   template <typename TSimType, template <typename> class TSimTraits> ImposedFieldReader<TSimType, TSimTraits>::ImposedFieldReader(typename TSimTraits<TSimType>::MagType &magB, typename TSimTraits<TSimType>::VelType &velV)
      : ImposedCodMagVelReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasVelICImposed, TSimTraits<TSimType>::HasVelOCImposed>("magVel", magB.oc().trunc())
   {
      this->mpMagB = &magB;

      this->mpVelV = &velV;
   }

   template <typename TSimType, template <typename> class TSimTraits> ImposedFieldReader<TSimType, TSimTraits>::ImposedFieldReader(typename TSimTraits<TSimType>::CodType &codC)
      : ImposedCodMagVelReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasVelICImposed, TSimTraits<TSimType>::HasVelOCImposed>("cod", codC.oc().trunc())
   {
      this->mpCodC = &codC;
   }

   template <typename TSimType, template <typename> class TSimTraits> ImposedFieldReader<TSimType, TSimTraits>::ImposedFieldReader(typename TSimTraits<TSimType>::MagType &magB)
      : ImposedCodMagVelReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasVelICImposed, TSimTraits<TSimType>::HasVelOCImposed>("mag", magB.oc().trunc())
   {
      this->mpMagB = &magB;
   }

   template <typename TSimType, template <typename> class TSimTraits> ImposedFieldReader<TSimType, TSimTraits>::ImposedFieldReader(typename TSimTraits<TSimType>::VelType &velV)
      : ImposedCodMagVelReader<TSimType, TSimTraits, TSimTraits<TSimType>::HasVelICImposed, TSimTraits<TSimType>::HasVelOCImposed>("vel", velV.oc().trunc())
   {
      this->mpVelV = &velV;
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
      this->readCod();

      // Read magnetic coefficients
      this->readMag();
      
      // Read velocity coefficients
      this->readVel();
   }

}

#endif // IMPOSEDFIELDREADER_HPP
