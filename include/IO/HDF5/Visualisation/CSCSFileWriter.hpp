/** \file CSCSFileWriter.hpp
 *  \brief Implementation of the CSCS visualisation file format writer
 */

#ifndef CSCSFILEWRITER_HPP
#define CSCSFILEWRITER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "IO/HDF5/Visualisation/CSCSFileWriterBase.hpp"
#include "Timestepping/TimestepParameters.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the CSCS visualisation file format writer
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class CSCSFileWriter: public CSCSFileWriterBase
   {
      public:
         /**
         * @brief Constructor with all three fields
         *
         * @param gridName HDF5 main group name
         * @param codC Codensity variable
         * @param magB Magnetic variable
         * @param velV Velocity variable
         * @param tsParams Timestep parameters
         */
         CSCSFileWriter(std::string gridName, const typename TSimTraits<TSimType>::CodType &codC, const typename TSimTraits<TSimType>::MagType &magB, const typename TSimTraits<TSimType>::VelType &velV, const TimestepParameters &tsParams);

         /**
         * @brief Constructor for thermal convection
         *
         * @param gridName HDF5 main group name
         * @param codC Codensity variable
         * @param velV Velocity variable
         * @param tsParams Timestep parameters
         */
         CSCSFileWriter(std::string gridName, const typename TSimTraits<TSimType>::CodType &codC, const typename TSimTraits<TSimType>::VelType &velV, const TimestepParameters &tsParams);

         /**
         * @brief Constructor for magneto-convection
         *
         * @param gridName HDF5 main group name
         * @param magB Magnetic variable
         * @param velV Velocity variable
         * @param tsParams Timestep parameters
         */
         CSCSFileWriter(std::string gridName, const typename TSimTraits<TSimType>::MagType &magB, const typename TSimTraits<TSimType>::VelType &velV, const TimestepParameters &tsParams);

         /**
         * @brief Constructor with only codensity
         *
         * @param gridName HDF5 main group name
         * @param codC Codensity variable
         * @param tsParams Timestep parameters
         */
         CSCSFileWriter(std::string gridName, const typename TSimTraits<TSimType>::CodType &codC, const TimestepParameters &tsParams);

         /**
         * @brief Constructor with only magnetic field
         *
         * @param gridName HDF5 main group name
         * @param magB Magnetic variable
         * @param tsParams Timestep parameters
         */
         CSCSFileWriter(std::string gridName, const typename TSimTraits<TSimType>::MagType &magB, const TimestepParameters &tsParams);

         /**
         * @brief Constructor for cases with all three fields
         *
         * @param gridName HDF5 main group name
         * @param velV Velocity variable
         * @param tsParams Timestep parameters
         */
         CSCSFileWriter(std::string gridName, const typename TSimTraits<TSimType>::VelType &velV, const TimestepParameters &tsParams);

         /**
         * @brief Destructor
         */
         virtual ~CSCSFileWriter() {};

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
          * @brief Reference to a TimestepParameters object to get access to time and timestep
          */
         const TimestepParameters&  mrTSParams;
   };

   template <typename TSimType, template <typename> class TSimTraits> CSCSFileWriter<TSimType, TSimTraits>::CSCSFileWriter(std::string gridName, const typename TSimTraits<TSimType>::CodType &codC, const typename TSimTraits<TSimType>::MagType &magB, const typename TSimTraits<TSimType>::VelType &velV, const TimestepParameters &tsParams)
      : CSCSFileWriterBase(gridName, codC.oc().trunc()), mpCodC(&codC), mpMagB(&magB), mpVelV(&velV), mrTSParams(tsParams)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> CSCSFileWriter<TSimType, TSimTraits>::CSCSFileWriter(std::string gridName, const typename TSimTraits<TSimType>::CodType &codC, const typename TSimTraits<TSimType>::VelType &velV, const TimestepParameters &tsParams)
      : CSCSFileWriterBase(gridName, codC.oc().trunc()), mpCodC(&codC), mpMagB(NULL), mpVelV(&velV), mrTSParams(tsParams)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> CSCSFileWriter<TSimType, TSimTraits>::CSCSFileWriter(std::string gridName, const typename TSimTraits<TSimType>::MagType &magB, const typename TSimTraits<TSimType>::VelType &velV, const TimestepParameters &tsParams)
      : CSCSFileWriterBase(gridName, magB.oc().trunc()), mpCodC(NULL), mpMagB(&magB), mpVelV(&velV), mrTSParams(tsParams)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> CSCSFileWriter<TSimType, TSimTraits>::CSCSFileWriter(std::string gridName, const typename TSimTraits<TSimType>::CodType &codC, const TimestepParameters &tsParams)
      : CSCSFileWriterBase(gridName, codC.oc().trunc()), mpCodC(&codC), mpMagB(NULL), mpVelV(NULL), mrTSParams(tsParams)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> CSCSFileWriter<TSimType, TSimTraits>::CSCSFileWriter(std::string gridName, const typename TSimTraits<TSimType>::MagType &magB, const TimestepParameters &tsParams)
      : CSCSFileWriterBase(gridName, magB.oc().trunc()), mpCodC(NULL), mpMagB(&magB), mpVelV(NULL), mrTSParams(tsParams)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> CSCSFileWriter<TSimType, TSimTraits>::CSCSFileWriter(std::string gridName, const typename TSimTraits<TSimType>::VelType &velV, const TimestepParameters &tsParams)
      : CSCSFileWriterBase(gridName, velV.oc().trunc()), mpCodC(NULL), mpMagB(NULL), mpVelV(&velV), mrTSParams(tsParams)
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void CSCSFileWriter<TSimType, TSimTraits>::write()
   {
      // Create file
      this->preWrite();

      // Create the header and version information
      this->createFileInfo();

      // Create the base group
      this->createBaseGroup();

      // Write real space grid to file
      this->writeGrid(this->mpTrunc->sim()->rad()->radGrid(), this->mpTrunc->sim()->hoz()->thGrid(), this->mpTrunc->sim()->hoz()->phGrid());

      // Write the codensity real space values
      if(this->mpCodC != NULL)
      {
         this->writeCodensity(this->mpCodC->oc().rtp().data());
      }

      // Write the magnetic coefficients
      if(this->mpMagB != NULL)
      {
         this->writeMagnetic(this->mpMagB->oc().rtp().r().data(), this->mpMagB->oc().rtp().theta().data(), this->mpMagB->oc().rtp().phi().data());
      }

      // Write the velocity coefficients
      if(this->mpVelV != NULL)
      {
         this->writeVelocity(this->mpVelV->oc().rtp().r().data(), this->mpVelV->oc().rtp().theta().data(), this->mpVelV->oc().rtp().phi().data());
      }

      // Close file
      this->postWrite();
   }

}

#endif // CSCSFILEWRITER_HPP
