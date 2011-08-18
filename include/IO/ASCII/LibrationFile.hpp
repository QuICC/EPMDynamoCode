/** \file LibrationFile.hpp
 *  \brief Implementation of a libration diagnostic file
 */

#ifndef LIBRATIONFILE_HPP
#define LIBRATIONFILE_HPP

// Configuration includes
//

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "IO/ASCII/LibrationFileDefs.hpp"
#include "IO/ASCII/ASCIIFieldWriter.hpp"
#include "IO/ASCII/ASCIIEWriter.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of a libration diagnostics file
    *
    * \param TFieldType Type of the field
    */
   template <typename TFieldType> class LibrationFile: public ASCIIFieldWriter<TFieldType, ASCIIEWriter>
   {
      public:
         /**
         * @brief Constructor
         *
         * @param var Field variable
         * @param name File name
         * @param tsParams Timestep parameters
         */
         LibrationFile(TFieldType& var, std::string name, const TimestepParameters &tsParams, const int format = 0);

         /**
         * @brief Destructor
         */
         virtual ~LibrationFile() {};

         /**
          * @brief Write energy to file
          */
         virtual void write();
         
      protected:

      private:
         /**
          * @brief Compute the poincare mode amplitude
          */
         void computePoincareMode();

         /**
          * @brief Setup the computation
          */
         void setup();

         /**
          * @brief Reference to a TimestepParameters object to get access to time and timestep
          */
         const TimestepParameters&  mrTSParams;

         /**
          * @brief Output format flag
          */
         const int mFormat;

         /**
          * @brief Storage for the Poincare mode amplitude
          */
         Array mAmplitude;

   };

   template <typename TFieldType> LibrationFile<TFieldType>::LibrationFile(TFieldType &var, std::string name, const TimestepParameters &tsParams, int format)
      : ASCIIFieldWriter<TFieldType, ASCIIEWriter>(var, name + LibrationFileDefs::BASENAME, LibrationFileDefs::EXTENSION, LibrationFileDefs::HEADER, LibrationFileDefs::VERSION), mrTSParams(tsParams), mFormat(format), mAmplitude(3)
   {
      // Setup the computation
      this->setup();
   }

   template <typename TFieldType> void LibrationFile<TFieldType>::write()
   {
      // Compute poincare mode
      this->computePoincareMode();

      if(this->doesIO())
      {
         //Do pre write operations
         this->preWrite();

         // Write energy
         this->mFile << this->mrTSParams.time() << "  " << this->mAmplitude.transpose();

         // Add newline at the end of ouput
         this->mFile << std::endl;

         //Do post write operations
         this->postWrite();
      }
   }

   template <typename TFieldType> void LibrationFile<TFieldType>::computePoincareMode()
   {
      // Get truncation information
      SmartTruncation pTrunc = this->mrVar.oc().trunc();

      this->mAmplitude.setConstant(0.0);

      // Get the stored ls
      ArrayI ls = this->trunc()->local()->spec()->lArray();
      ArrayI ms;
      int l_ = -1;
      int m_ = -1;

      // Setup the x axis solid body rotation
      ArrayZ xAxis = ArrayZ::Zero(pTrunc->sim()->rad()->nN());
      xAxis(0).real() = -0.443113;
      EPMFloat xProj = 0.0;

      // Setup the y axis solid body rotation
      ArrayZ yAxis = ArrayZ::Zero(pTrunc->sim()->rad()->nN());
      yAxis(0).imag() = 0.443113;
      EPMFloat yProj = 0.0;

      // Loop over all stored l's
      for(int l = 0; l < ls.size(); l++)
      {
         l_ = ls(l);
         if(l_ == 1)
         {
            ms = this->trunc()->local()->spec()->mArray(l);
            for(int m =0; m < ms.size(); m++)
            {
               m_ = ms(m);
               if(m_ == 1)
               {
                  EPMFloat shFactor = 4.0*4.0 * MathConstants::PI / static_cast<EPMFloat>(2*l_+1);
                  EPMFloat lfactor = static_cast<EPMFloat>(l_*(l_+1));
                  int nN = pTrunc->sim()->rad()->nN();

                  for(int n = 0; n < nN; ++n)
                  {
                     for(int k = 0; k < nN; ++k)
                     {
                        xProj += this->mTransform.radBasis().at(l).eWeights()(k,n) * (this->velV().oc().perturbation().tor().lshell(l)(n,m).real()*xAxis(k).real() + this->velV().oc().perturbation().tor().lshell(l)(n,m).imag()*xAxis(k).imag());
                        yProj += this->mTransform.radBasis().at(l).eWeights()(k,n) * (this->velV().oc().perturbation().tor().lshell(l)(n,m).real()*yAxis(k).real() + this->velV().oc().perturbation().tor().lshell(l)(n,m).imag()*yAxis(k).imag());
                     }
                  }

                  xProj *= shFactor*lfactor;
                  yProj *= shFactor*lfactor;

                  this->mAmplitude(0) += xProj;
                  this->mAmplitude(1) += yProj;
               }
            }
         }
      }

      #ifdef EPMDYNAMO_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->mAmplitude.data(), this->mAmplitude.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif // EPMDYNAMO_MPI
   }

}

#endif // LIBRATIONFILE_HPP
