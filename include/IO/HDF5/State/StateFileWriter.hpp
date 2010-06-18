/** \file StateFileWriter.hpp
 *  \brief Implementation of the HDF5 state file writer
 */

#ifndef STATEFILEWRITER_HPP
#define STATEFILEWRITER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "IO/HDF5/State/StateFileWriterBase.hpp"
#include "Timestepping/TimestepParameters.hpp"
#include "Equations/Parameters/EquationParameters.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the HDF5 state file writer
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits of the implementation
    */
   template <typename TSimType, template <typename> class TSimTraits> class StateFileWriter: public StateFileWriterBase
   {
      public:
         /**
         * @brief Constructor with all three fields
         *
         * @param codC Codensity variable
         * @param magB Magnetic variable
         * @param velV Velocity variable
         * @param eqParams Equation parameters
         * @param tsParams Timestep parameters
         */
         StateFileWriter(const typename TSimTraits<TSimType>::CodType &codC, const typename TSimTraits<TSimType>::MagType &magB, const typename TSimTraits<TSimType>::VelType &velV, const EquationParameters &eqParams, const TimestepParameters &tsParams);

         /**
         * @brief Constructor for thermal convection
         *
         * @param codC Codensity variable
         * @param velV Velocity variable
         * @param eqParams Equation parameters
         * @param tsParams Timestep parameters
         */
         StateFileWriter(const typename TSimTraits<TSimType>::CodType &codC, const typename TSimTraits<TSimType>::VelType &velV, const EquationParameters &eqParams, const TimestepParameters &tsParams);

         /**
         * @brief Constructor for magneto-convection
         *
         * @param magB Magnetic variable
         * @param velV Velocity variable
         * @param eqParams Equation parameters
         * @param tsParams Timestep parameters
         */
         StateFileWriter(const typename TSimTraits<TSimType>::MagType &magB, const typename TSimTraits<TSimType>::VelType &velV, const EquationParameters &eqParams, const TimestepParameters &tsParams);

         /**
         * @brief Constructor with only codensity
         *
         * @param codC Codensity variable
         * @param eqParams Equation parameters
         * @param tsParams Timestep parameters
         */
         StateFileWriter(const typename TSimTraits<TSimType>::CodType &codC, const EquationParameters &eqParams, const TimestepParameters &tsParams);

         /**
         * @brief Constructor with only magnetic field
         *
         * @param magB Magnetic variable
         * @param eqParams Equation parameters
         * @param tsParams Timestep parameters
         */
         StateFileWriter(const typename TSimTraits<TSimType>::MagType &magB, const EquationParameters &eqParams, const TimestepParameters &tsParams);

         /**
         * @brief Constructor for cases with all three fields
         *
         * @param velV Velocity variable
         * @param eqParams Equation parameters
         * @param tsParams Timestep parameters
         */
         StateFileWriter(const typename TSimTraits<TSimType>::VelType &velV, const EquationParameters &eqParams, const TimestepParameters &tsParams);

         /**
         * @brief Destructor
         */
         virtual ~StateFileWriter() {};

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

         /**
          * @brief Store physical parameters in a array for faster access
          */
         Array   mPhys;

         /**
          * @brief Reference to a TimestepParameters object to get access to time and timestep
          */
         const TimestepParameters&  mrTSParams;

   };

   template <typename TSimType, template <typename> class TSimTraits> StateFileWriter<TSimType, TSimTraits>::StateFileWriter(const typename TSimTraits<TSimType>::CodType &codC, const typename TSimTraits<TSimType>::MagType &magB, const typename TSimTraits<TSimType>::VelType &velV, const EquationParameters &eqParams, const TimestepParameters &tsParams)
      : StateFileWriterBase(codC.trunc()), mpCodC(&codC), mpMagB(&magB), mpVelV(&velV), mTrunc(4), mPhys(4), mrTSParams(tsParams)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();

      this->mPhys(0) = eqParams.E();
      this->mPhys(1) = eqParams.q();
      this->mPhys(2) = eqParams.Ra();
      this->mPhys(3) = eqParams.Ro();
   }

   template <typename TSimType, template <typename> class TSimTraits> StateFileWriter<TSimType, TSimTraits>::StateFileWriter(const typename TSimTraits<TSimType>::CodType &codC, const typename TSimTraits<TSimType>::VelType &velV, const EquationParameters &eqParams, const TimestepParameters &tsParams)
      : StateFileWriterBase(codC.trunc()), mpCodC(&codC), mpMagB(NULL), mpVelV(&velV), mTrunc(4), mPhys(4), mrTSParams(tsParams)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();

      this->mPhys(0) = eqParams.E();
      this->mPhys(1) = eqParams.q();
      this->mPhys(2) = eqParams.Ra();
      this->mPhys(3) = eqParams.Ro();
   }

   template <typename TSimType, template <typename> class TSimTraits> StateFileWriter<TSimType, TSimTraits>::StateFileWriter(const typename TSimTraits<TSimType>::MagType &magB, const typename TSimTraits<TSimType>::VelType &velV, const EquationParameters &eqParams, const TimestepParameters &tsParams)
      : StateFileWriterBase(magB.trunc()), mpCodC(NULL), mpMagB(&magB), mpVelV(&velV), mTrunc(4), mPhys(4), mrTSParams(tsParams)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();

      this->mPhys(0) = eqParams.E();
      this->mPhys(1) = eqParams.q();
      this->mPhys(2) = eqParams.Ra();
      this->mPhys(3) = eqParams.Ro();
   }

   template <typename TSimType, template <typename> class TSimTraits> StateFileWriter<TSimType, TSimTraits>::StateFileWriter(const typename TSimTraits<TSimType>::CodType &codC, const EquationParameters &eqParams, const TimestepParameters &tsParams)
      : StateFileWriterBase(codC.trunc()), mpCodC(&codC), mpMagB(NULL), mpVelV(NULL), mTrunc(4), mPhys(4), mrTSParams(tsParams)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();

      this->mPhys(0) = eqParams.E();
      this->mPhys(1) = eqParams.q();
      this->mPhys(2) = eqParams.Ra();
      this->mPhys(3) = eqParams.Ro();
   }

   template <typename TSimType, template <typename> class TSimTraits> StateFileWriter<TSimType, TSimTraits>::StateFileWriter(const typename TSimTraits<TSimType>::MagType &magB, const EquationParameters &eqParams, const TimestepParameters &tsParams)
      : StateFileWriterBase(magB.trunc()), mpCodC(NULL), mpMagB(&magB), mpVelV(NULL), mTrunc(4), mPhys(4), mrTSParams(tsParams)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();

      this->mPhys(0) = eqParams.E();
      this->mPhys(1) = eqParams.q();
      this->mPhys(2) = eqParams.Ra();
      this->mPhys(3) = eqParams.Ro();
   }

   template <typename TSimType, template <typename> class TSimTraits> StateFileWriter<TSimType, TSimTraits>::StateFileWriter(const typename TSimTraits<TSimType>::VelType &velV, const EquationParameters &eqParams, const TimestepParameters &tsParams)
      : StateFileWriterBase(velV.trunc()), mpCodC(NULL), mpMagB(NULL), mpVelV(&velV), mTrunc(4), mPhys(4), mrTSParams(tsParams)
   {
      this->mTrunc(0) = this->mpTrunc->sim()->rad()->maxN();
      this->mTrunc(1) = this->mpTrunc->sim()->hoz()->maxL();
      this->mTrunc(2) = this->mpTrunc->sim()->hoz()->maxM();
      this->mTrunc(3) = this->mpTrunc->sim()->hoz()->mp();

      this->mPhys(0) = eqParams.E();
      this->mPhys(1) = eqParams.q();
      this->mPhys(2) = eqParams.Ra();
      this->mPhys(3) = eqParams.Ro();
   }

   template <typename TSimType, template <typename> class TSimTraits> void StateFileWriter<TSimType, TSimTraits>::write()
   {
      // Create file
      this->preWrite();

      // Create the header and version information
      this->createFileInfo();

      // Write the Physical parameters
      this->writePhysical(this->mPhys(0), this->mPhys(1), this->mPhys(2), this->mPhys(3));

      // Write the truncation information
      this->writeTruncation(this->mTrunc(0), this->mTrunc(1), this->mTrunc(2), this->mTrunc(3));

      // Write the run information
      this->writeRun(this->mrTSParams.time(), this->mrTSParams.dt());

      // Write the codensity coefficients
      if(TSimTraits<TSimType>::NeedCodensity)
      {
         this->writeCodensity(this->mpCodC->oc().perturbation().data());
      }

      // Write the magnetic coefficients
      if(TSimTraits<TSimType>::NeedMagnetic)
      {
         this->writeMagnetic(this->mpMagB->oc().perturbation().tor().data(), this->mpMagB->oc().perturbation().pol().data());
      }

      // Write the velocity coefficients
      if(TSimTraits<TSimType>::NeedVelocity)
      {
         this->writeVelocity(this->mpVelV->oc().perturbation().tor().data(), this->mpVelV->oc().perturbation().pol().data());
      }

      // Close file
      this->postWrite();
   }

}

#endif // STATEFILEWRITER_HPP
