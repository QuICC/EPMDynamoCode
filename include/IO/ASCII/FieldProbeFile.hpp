/** \file FieldProbeFile.hpp
 *  \brief Implementation of an field probe file
 */

#ifndef FIELDPROBEFILE_HPP
#define FIELDPROBEFILE_HPP

// Configuration includes
//

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "IO/ASCII/FieldProbeFileDefs.hpp"
#include "IO/ASCII/ASCIIFieldWriter.hpp"
#include "IO/ASCII/ASCIIEWriter.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of ASCIIFile for a field probe file
    *
    * \param TFieldType Type of the field
    */
   template <typename TFieldType> class FieldProbeFile: public ASCIIFieldWriter<TFieldType, ASCIIEWriter>
   {
      public:
         /**
          * @brief Enum of the possible probes
          */
         enum FieldProbe {RADIALPROBE, THETAPROBE, PHIPROBE, FULLPROBE};

         /**
          * @brief Constructor
          *
          * @param var Field variable
          * @param name File name
          * @param posRatios Ratios of full grid defining position
          * @param tsParams Timestep parameters
          */
         FieldProbeFile(TFieldType& var, std::string name, const TimestepParameters &tsParams, const Array& posRatios, FieldProbe probe = FULLPROBE);

         /**
          * @brief Destructor
          */
         virtual ~FieldProbeFile() {};

         /**
          * @brief Write energy to file
          */
         virtual void write();
         
      protected:

      private:
         /**
          * @brief Initialise the probe
          *
          * @param posRatios Ratios of full grid defining position
          */
         void initProbe(const Array& posRatios);

         /**
          * @brief Initialise the probe
          */
         void updateProbe();

         /**
          * @brief Reference to a TimestepParameters object to get access to time and timestep
          */
         const TimestepParameters&  mrTSParams;

         /**
          * @brief Output probe flag
          */
         FieldProbe mProbe;

         /**
          * @brief Position of the probe (simulation wide)
          */
         ArrayI mSimPos;

         /**
          * @brief Position of the probe (local)
          */
         ArrayI mPos;

         /**
          * @brief Storage for the data
          */
         Array mData;
   };

   template <typename TFieldType> FieldProbeFile<TFieldType>::FieldProbeFile(TFieldType &var, std::string name, const TimestepParameters &tsParams, const Array& posRatios, FieldProbe probe)
      : ASCIIFieldWriter<TFieldType, ASCIIEWriter>(var, name + FieldProbeFileDefs::BASENAME, FieldProbeFileDefs::EXTENSION, FieldProbeFileDefs::HEADER, FieldProbeFileDefs::VERSION), mrTSParams(tsParams), mProbe(probe), mSimPos(3), mPos(3), mData(1)
   {
      this->initProbe(posRatios);
   }

   template <typename TFieldType> void FieldProbeFile<TFieldType>::initProbe(const Array& posRatios)
   {
      // Resize data storage if full probe is required
      if(this->mProbe == FULLPROBE)
      {
         this->mData.resize(3);
      }

      SmartTruncation pTrunc = this->mrVar.oc().trunc();

      // Set radial position
      this->mSimPos(0) = pTrunc->sim()->rad()->nR()*posRatios(0);
      // Set theta position
      this->mSimPos(1) = pTrunc->sim()->hoz()->nTh()*posRatios(1);
      // Set phi position
      this->mSimPos(2) = pTrunc->sim()->hoz()->nPh()*posRatios(2);

      EPMFloat radius = pTrunc->sim()->rad()->radGrid()(this->mSimPos(0));
      EPMFloat cosTheta = pTrunc->sim()->hoz()->cosTheta()(this->mSimPos(1));

      this->mPos.setConstant(-1);

      for(int n = 0; n < pTrunc->local()->rtp()->nR(); n++)
      {
         if(radius == pTrunc->local()->rtp()->radGrid(n))
         {
            this->mPos(0) = n;

            for(int th = 0; th < pTrunc->local()->rtp()->nTh(n); th++)
            {
               if(cosTheta == pTrunc->local()->rtp()->cTh(th, n))
               {
                  this->mPos(1) = th;
                  this->mPos(2) = this->mSimPos(2);
                  break;
               }
            }

            break;
         }
      }

      if(this->mPos(2) == -1)
      {
         this->mPos.resize(0);
      }
   }

   template <typename TFieldType> void FieldProbeFile<TFieldType>::updateProbe()
   {
      if(this->mPos.size() == 3)
      {
         if(this->mProbe == RADIALPROBE)
         {
            this->mData(0) = this->mrVar.oc().rtp().r().shell(this->mPos(0))(this->mPos(2), this->mPos(1));
         }
         else if(this->mProbe == THETAPROBE)
         {
            this->mData(0) = this->mrVar.oc().rtp().theta().shell(this->mPos(0))(this->mPos(2), this->mPos(1));
         }
         // If format = 2 output the energy harmonic order spectrum evolution
         else if(this->mProbe == PHIPROBE)
         {
            this->mData(0) = this->mrVar.oc().rtp().phi().shell(this->mPos(0))(this->mPos(2), this->mPos(1));
         }
         // If format = 3 output the symmetric/anti-symmetric energy decomposition evolution
         else if(this->mProbe == FULLPROBE)
         {
            this->mData(0) = this->mrVar.oc().rtp().r().shell(this->mPos(0))(this->mPos(2), this->mPos(1));
            this->mData(1) = this->mrVar.oc().rtp().theta().shell(this->mPos(0))(this->mPos(2), this->mPos(1));
            this->mData(2) = this->mrVar.oc().rtp().phi().shell(this->mPos(0))(this->mPos(2), this->mPos(1));
         }
      } else
      {
         this->mData.setConstant(0.0);
      }

      // Get the "global" probe from MPI code
      #ifdef EPMDYNAMO_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->mData.data(), this->mData.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif // EPMDYNAMO_MPI
   }


   template <typename TFieldType> void FieldProbeFile<TFieldType>::write()
   {
      // Update the data over parallelisation
      this->updateProbe();

      if(this->doesIO())
      {
         //Do pre write operations
         this->preWrite();

         SmartTruncation pTrunc = this->mrVar.oc().trunc();

         EPMFloat radius = pTrunc->sim()->rad()->radGrid()(this->mSimPos(0));
         EPMFloat theta = pTrunc->sim()->hoz()->thGrid()(this->mSimPos(1));
         EPMFloat phi = pTrunc->sim()->hoz()->phGrid()(this->mSimPos(2));

         // Write energy
         this->mFile << std::setprecision(16) << this->mrTSParams.time() << "  " << radius << "  " << theta << "  " << phi << "  " << this->mData.transpose();

         // Add newline at the end of ouput
         this->mFile << std::endl;

         //Do post write operations
         this->postWrite();
      }
   }

}

#endif // FIELDPROBEFILE_HPP
