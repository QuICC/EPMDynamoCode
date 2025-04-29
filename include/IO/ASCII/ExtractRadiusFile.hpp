/** \file ExtractRadiusFile.hpp
 *  \brief Implementation of a libration diagnostic file
 */

#ifndef EXTRACTRADIUSFILE_HPP
#define EXTRACTRADIUSFILE_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"

// System includes
//
#include <string>
#include <iomanip>

// External includes
//

// Project includes
//
#include "General/MathConstants.hpp"
#include "IO/ASCII/ExtractRadiusFileDefs.hpp"
#include "IO/ASCII/ASCIIFieldWriter.hpp"
#include "IO/ASCII/ASCIIEWriter.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of a libration diagnostics file
    *
    * \param TFieldType Type of the field
    */
   template <typename TFieldType> class ExtractRadiusFile: public ASCIIFieldWriter<TFieldType, ASCIIEWriter>
   {
      public:
         /// Typedef for the Spectral radial Transform data type
         typedef typename SimulationConfig::NumericalScheme::RadialTransformType::PolynomialBasis RadialBasisType;

         /**
         * @brief Constructor
         *
         * @param var Field variable
         * @param name File name
         */
         ExtractRadiusFile(TFieldType& var, const int theta, const int phi, std::string name, const int format = 0);

         /**
         * @brief Destructor
         */
         virtual ~ExtractRadiusFile() {};

         /**
          * @brief Write energy to file
          */
         virtual void write();
         
      protected:

      private:
         /**
          * @brief Gather the data along the radius
          */
         void gatherRadius();
          
         /**
          * @brief Output format flag
          */
         const int mFormat;

         /**
          * @brief Theta point at which to take radius
          */
         int mTheta;

         /**
          * @brief Phi point at which to take radius
          */
         int mPhi;

         /**
          * @brief Storage for the Poincare mode amplitude
          */
         Matrix mRadius;

   };

   template <typename TFieldType> ExtractRadiusFile<TFieldType>::ExtractRadiusFile(TFieldType &var, const int theta, const int phi, std::string name, int format)
      : ASCIIFieldWriter<TFieldType, ASCIIEWriter>(var, name + ExtractRadiusFileDefs::BASENAME, ExtractRadiusFileDefs::EXTENSION, ExtractRadiusFileDefs::HEADER, ExtractRadiusFileDefs::VERSION), mFormat(format), mTheta(theta), mPhi(phi), mRadius(var.oc().trunc()->sim()->rad()->nR(),1)
   {
      if(this->mPhi < 0)
      {
         this->mRadius.resize(this->mrVar.oc().trunc()->sim()->rad()->nR(), this->mrVar.oc().trunc()->sim()->hoz()->nPh());
         this->mPhi = 0;
      }
   }

   template <typename TFieldType> void ExtractRadiusFile<TFieldType>::write()
   {
      // Gather the data
      this->gatherRadius();

      if(this->doesIO())
      {
         //Do pre write operations
         this->preWrite();

         this->mFile << "# Selected latitude: " << (this->mrVar.oc().trunc()->sim()->hoz()->thGrid()(this->mTheta)/MathConstants::PI)*180. << std::endl;
         if(this->mRadius.cols() > 1)
         {
            this->mFile << "# All longitudes have been extracted" << std::endl;
            this->mFile << "# r min max phi_0 ... phi_n" << std::endl;
         } else
         {
            this->mFile << "# Selected longitude: " << (this->mrVar.oc().trunc()->sim()->hoz()->phGrid()(this->mPhi)/MathConstants::PI)*180. << std::endl;
         }

         // Write radius values
         for(int r = 0; r < this->mRadius.rows(); r++)
         {
            this->mFile << this->mrVar.oc().trunc()->sim()->rad()->radGrid()(r) << "  ";

            // Compute the envelope
            if(this->mRadius.cols() > 1)
            {
               this->mFile << this->mRadius.row(r).minCoeff() << "  " << this->mRadius.row(r).maxCoeff() << "  ";
            }

            this->mFile << this->mRadius.row(r) << std::endl;
         }

         // Add newline at the end of ouput
         this->mFile << std::endl;

         //Do post write operations
         this->postWrite();
      }
   }

   template <typename TFieldType> void ExtractRadiusFile<TFieldType>::gatherRadius()
   {
      // Get the simulation wide first radial point
      int r0 = this->mrVar.oc().trunc()->local()->rtp()->r0();

      this->mRadius.setConstant(0.0);
      int th;
      if(this->mFormat == 0)
      {
         for(int r = 0; r < this->mrVar.oc().trunc()->local()->rtp()->nR(); r++)
         {
            if(this->mTheta > this->mrVar.oc().trunc()->local()->rtp()->th0(r) - 1 && this->mTheta < this->mrVar.oc().trunc()->local()->rtp()->th0(r) + this->mrVar.oc().trunc()->local()->rtp()->nTh(r))
            {
               th = this->mTheta - this->mrVar.oc().trunc()->local()->rtp()->th0(r);
               this->mRadius.row(r0 + r) = this->mrVar.oc().rtp().r().shell(r).col(th).segment(this->mPhi, this->mRadius.cols()).transpose();
            }
         }
      } else if(this->mFormat == 1)
      {
         for(int r = 0; r < this->mrVar.oc().trunc()->local()->rtp()->nR(); r++)
         {
            if(this->mTheta > this->mrVar.oc().trunc()->local()->rtp()->th0(r) - 1 && this->mTheta < this->mrVar.oc().trunc()->local()->rtp()->th0(r) + this->mrVar.oc().trunc()->local()->rtp()->nTh(r))
            {
               th = this->mTheta - this->mrVar.oc().trunc()->local()->rtp()->th0(r);
               this->mRadius.row(r0 + r) = this->mrVar.oc().rtp().theta().shell(r).col(th).segment(this->mPhi, this->mRadius.cols()).transpose();
            }
         }
      } else if(this->mFormat == 2)
      {
         for(int r = 0; r < this->mrVar.oc().trunc()->local()->rtp()->nR(); r++)
         {
            if(this->mTheta > this->mrVar.oc().trunc()->local()->rtp()->th0(r) - 1 && this->mTheta < this->mrVar.oc().trunc()->local()->rtp()->th0(r) + this->mrVar.oc().trunc()->local()->rtp()->nTh(r))
            {
               th = this->mTheta - this->mrVar.oc().trunc()->local()->rtp()->th0(r);
               this->mRadius.row(r0 + r) = this->mrVar.oc().rtp().phi().shell(r).col(th).segment(this->mPhi, this->mRadius.cols()).transpose();
            }
         }
      }

      // Get the "global" radius from MPI code
      #ifdef EPMDYNAMO_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->mRadius.data(), this->mRadius.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif // EPMDYNAMO_MPI
   }

}

#endif // EXTRACTRADIUSFILE_HPP
