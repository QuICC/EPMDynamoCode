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
#include "IO/HDF5/Visualisation/CSCSFileDefs.hpp"
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
          * @brief Write basic data to file
          */
         virtual void write();

         /**
          * @brief Write general information to file
          */
         void writeStart();

         /**
          * @brief Finish writting to file
          */
         void writeEnd();

         /**
          * @brief Write visualisation data to file
          *
          * @param filter Filter tag name
          */
         template <typename TVisTraits> void writeVisualisation(const std::string& filter);
         
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

      // Close file
      this->postWrite();
   }

   template <typename TSimType, template <typename> class TSimTraits> void CSCSFileWriter<TSimType, TSimTraits>::writeStart()
   {
      // Create file
      this->preWrite();

      // Create the header and version information
      this->createFileInfo();

      // Create the base group
      this->createBaseGroup();

      // Write real space grid to file
      this->writeGrid(this->mpTrunc->sim()->rad()->radGrid(), this->mpTrunc->sim()->hoz()->thGrid(), this->mpTrunc->sim()->hoz()->phGrid());
   }

   template <typename TSimType, template <typename> class TSimTraits> void CSCSFileWriter<TSimType, TSimTraits>::writeEnd()
   {
      // Close file
      this->postWrite();
   }

   template <typename TSimType, template <typename> class TSimTraits> template <typename TVisTraits> void CSCSFileWriter<TSimType, TSimTraits>::writeVisualisation(const std::string&  filter = "")
   {
      // Write the codensity values
      if(this->mpCodC != NULL && TVisTraits::VisCodensity)
      {
         this->writeScalarField(CSCSFileDefs::CODENSITYTAG+filter, this->mpCodC->oc().rtp().data());
      }

      // Write the codensity gradient values
      if(this->mpCodC != NULL && TVisTraits::VisCodensityGrad)
      {
         this->writeVectorField(CSCSFileDefs::CODENSITYGRADTAG+filter, this->mpCodC->oc().grad().r().data(), this->mpCodC->oc().grad().theta().data(), this->mpCodC->oc().grad().phi().data());
      }

      // Write the magnetic field values
      if(this->mpMagB != NULL && TVisTraits::VisMagnetic)
      {
         this->writeVectorField(CSCSFileDefs::MAGNETICTAG+filter, this->mpMagB->oc().rtp().r().data(), this->mpMagB->oc().rtp().theta().data(), this->mpMagB->oc().rtp().phi().data());
      }

      // Write the Lorentz force values
      if(this->mpMagB != NULL && TVisTraits::VisLorentz)
      {
         RTPField tmp(this->mpVelV->oc().rtp().trunc());
         this->mpMagB->oc().rtp().template cross<0>(tmp, this->mpMagB->oc().curl());
         this->writeVectorField(CSCSFileDefs::LORENTZTAG+filter, tmp.r().data(), tmp.theta().data(), tmp.phi().data());
      }

      // Write the velocity field values
      if(this->mpVelV != NULL && TVisTraits::VisVelocity)
      {
         this->writeVectorField(CSCSFileDefs::VELOCITYTAG+filter, this->mpVelV->oc().rtp().r().data(), this->mpVelV->oc().rtp().theta().data(), this->mpVelV->oc().rtp().phi().data());
      }

      // Write the vorticity values
      if(this->mpVelV != NULL && TVisTraits::VisVorticity)
      {
         this->writeVectorField(CSCSFileDefs::VORTICITYTAG+filter, this->mpVelV->oc().curl().r().data(), this->mpVelV->oc().curl().theta().data(), this->mpVelV->oc().curl().phi().data());
      }

      // Write the vorticity values
      if(this->mpVelV != NULL && TVisTraits::VisHelicity)
      {
         RTPScalar tmp(this->mpVelV->oc().rtp().trunc());
         this->mpVelV->oc().rtp().template dot<0>(tmp, this->mpVelV->oc().curl(), 1.0);
         this->writeScalarField(CSCSFileDefs::HELICITYTAG+filter, tmp.data());
      }
   }

}

#endif // CSCSFILEWRITER_HPP
