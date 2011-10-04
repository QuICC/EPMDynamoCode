/** \file StateFileReader.hpp
 *  \brief Implementation of HDF5 state file reader
 */

#ifndef STATEFILEREADER_HPP
#define STATEFILEREADER_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "IO/HDF5/State/StateFileDefs.hpp"
#include "IO/HDF5/State/StateFileReaderBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of HDF5 state file reader
    *
    * \tparam TSimTraits Traits of the implementation
    */
   template <typename TSimTraits> class StateFileReader: public StateFileReaderBase
   {
      public:
         /**
         * @brief Constructor with all three fields
         *
         * @param codC Codensity variable
         * @param magB Magnetic variable
         * @param velV Velocity variable
         * @param name file name
         */
         StateFileReader(typename TSimTraits::CodType& codC, typename TSimTraits::MagType& magB, typename TSimTraits::VelType& velV, std::string name);

         /**
         * @brief Constructor for thermal convection
         *
         * @param codC Codensity variable
         * @param velV Velocity variable
         * @param name file name
         */
         StateFileReader(typename TSimTraits::CodType& codC, typename TSimTraits::VelType& velV, std::string name);

         /**
         * @brief Constructor for magneto-convection
         *
         * @param magB Magnetic variable
         * @param velV Velocity variable
         * @param name file name
         */
         StateFileReader(typename TSimTraits::MagType& magB, typename TSimTraits::VelType& velV, std::string name);

         /**
         * @brief Constructor with only codensity
         *
         * @param codC Codensity variable
         * @param name file name
         */
         StateFileReader(typename TSimTraits::CodType& codC, std::string name);

         /**
         * @brief Constructor with only magnetic field
         *
         * @param magB Magnetic variable
         * @param name file name
         */
         StateFileReader(typename TSimTraits::MagType& magB, std::string name);

         /**
         * @brief Constructor for cases with all three fields
         *
         * @param velV Velocity variable
         * @param name file name
         */
         StateFileReader(typename TSimTraits::VelType& velV, std::string name);

         /**
         * @brief Destructor
         */
         virtual ~StateFileReader() {};

         /**
          * @brief Read State from file
          */
         virtual void read();

         /**
          * @brief Read resolution and data information from state file
          */
         void readSetup();

         /**
          * @brief Read Partial data from state file
          */
         template <typename TFilter> void readPartial();
         
      protected:
         /**
          * @brief Pointer to the condensity variable
          */
         typename TSimTraits::CodType*  mpCodC;

         /**
          * @brief Pointer to the magnetic variable
          */
         typename TSimTraits::MagType*  mpMagB;

         /**
          * @brief Pointer to the velocity variable
          */
         typename TSimTraits::VelType*  mpVelV;

      private:
   };

   template <typename TSimTraits> StateFileReader<TSimTraits>::StateFileReader(typename TSimTraits::CodType &codC, typename TSimTraits::MagType &magB, typename TSimTraits::VelType &velV, std::string name)
      : StateFileReaderBase(name, codC.oc().trunc()), mpCodC(&codC), mpMagB(&magB), mpVelV(&velV)
   {
      // Make sur the extra terms are zero
      this->mpCodC->initialiseZeros();

      // Make sur the extra terms are zero
      this->mpMagB->initialiseZeros();

      // Make sur the extra terms are zero
      this->mpVelV->initialiseZeros();
   }

   template <typename TSimTraits> StateFileReader<TSimTraits>::StateFileReader(typename TSimTraits::CodType &codC, typename TSimTraits::VelType &velV, std::string name)
      : StateFileReaderBase(name, codC.oc().trunc()), mpCodC(&codC), mpMagB(NULL), mpVelV(&velV)
   {
      // Make sur the extra terms are zero
      this->mpCodC->initialiseZeros();

      // Make sur the extra terms are zero
      this->mpVelV->initialiseZeros();
   }

   template <typename TSimTraits> StateFileReader<TSimTraits>::StateFileReader(typename TSimTraits::MagType &magB, typename TSimTraits::VelType &velV, std::string name)
      : StateFileReaderBase(name, magB.oc().trunc()), mpCodC(NULL), mpMagB(&magB), mpVelV(&velV)
   {
      // Make sur the extra terms are zero
      this->mpMagB->initialiseZeros();

      // Make sur the extra terms are zero
      this->mpVelV->initialiseZeros();
   }

   template <typename TSimTraits> StateFileReader<TSimTraits>::StateFileReader(typename TSimTraits::CodType &codC, std::string name)
      : StateFileReaderBase(name, codC.oc().trunc()), mpCodC(&codC), mpMagB(NULL), mpVelV(NULL)
   {
      // Make sur the extra terms are zero
      this->mpCodC->initialiseZeros();
   }

   template <typename TSimTraits> StateFileReader<TSimTraits>::StateFileReader(typename TSimTraits::MagType &magB, std::string name)
      : StateFileReaderBase(name, magB.oc().trunc()), mpCodC(NULL), mpMagB(&magB), mpVelV(NULL)
   {
      // Make sur the extra terms are zero
      this->mpMagB->initialiseZeros();
   }

   template <typename TSimTraits> StateFileReader<TSimTraits>::StateFileReader(typename TSimTraits::VelType &velV, std::string name)
      : StateFileReaderBase(name, velV.oc().trunc()), mpCodC(NULL), mpMagB(NULL), mpVelV(&velV)
   {
      // Make sur the extra terms are zero
      this->mpVelV->initialiseZeros();
   }

   template <typename TSimTraits> void StateFileReader<TSimTraits>::read()
   {
      // THIS IS NOT GENERAL ENOUGH ANYMORE, BUT ALSO NOT USEFUL IN THE CURRENT STATE
      // Read the Physical parameters
      // EPMFloat E, q, Ra, Ro;
      // this->readPhysical(E, q, Ra, Ro);

      // Read the truncation information
      this->readTruncation();

      // Check file compatibility with data truncation
      this->checkTruncation();

      // Set Read arguments
      this->setReadArguments();

      // Read the run information
      this->readRun();

      // Read codensity coefficients
      if(TSimTraits::NeedCodensity)
      {
         this->readScalarField(StateFileDefs::CODENSITYTAG, this->mpCodC->rOc().rPerturbation().data());

         this->mpCodC->rOc().rPerturbation().checkInitialisation();
      }

      // Read magnetic coefficients
      if(TSimTraits::NeedMagnetic)
      {
         this->readTorPolField(StateFileDefs::MAGNETICTAG, this->mpMagB->rOc().rPerturbation().rTor().data(), this->mpMagB->rOc().rPerturbation().rPol().data());

         this->mpMagB->rOc().rPerturbation().checkInitialisation();
      }
      
      // Read velocity coefficients
      if(TSimTraits::NeedVelocity)
      {
         this->readTorPolField(StateFileDefs::VELOCITYTAG, this->mpVelV->rOc().rPerturbation().rTor().data(), this->mpVelV->rOc().rPerturbation().rPol().data());

         this->mpVelV->rOc().rPerturbation().checkInitialisation();
      }
   }

   template <typename TSimTraits> void StateFileReader<TSimTraits>::readSetup()
   {
      // THIS IS NOT GENERAL ENOUGH ANYMORE, BUT ALSO NOT USEFUL IN THE CURRENT STATE
      // Read the Physical parameters
      // EPMFloat E, q, Ra, Ro;
      // this->readPhysical(E, q, Ra, Ro);

      // Read the truncation information
      this->readTruncation();

      // Check file compatibility with data truncation
      this->checkTruncation();

      // Set Read arguments
      this->setReadArguments();

      // Read the run information
      this->readRun();
   }

   template <typename TSimTraits> template <typename TFilter> void StateFileReader<TSimTraits>::readPartial()
   {
      // Read requested codensity coefficients
      if(TSimTraits::NeedCodensity && TFilter::ReadCodensity == StateFileDefs::FullField)
      {
         this->readScalarField(StateFileDefs::CODENSITYTAG, this->mpCodC->rOc().rPerturbation().data());
      }

      // Read requested magnetic coefficients
      if(TSimTraits::NeedMagnetic && TFilter::ReadMagnetic == StateFileDefs::FullField)
      {
         this->readTorPolField(StateFileDefs::MAGNETICTAG, this->mpMagB->rOc().rPerturbation().rTor().data(), this->mpMagB->rOc().rPerturbation().rPol().data());
      } else if(TSimTraits::NeedMagnetic && TFilter::ReadMagnetic == StateFileDefs::ToroidalOnly)
      {
         // Read Toroidal component
         this->readTorField(StateFileDefs::MAGNETICTAG, this->mpMagB->rOc().rPerturbation().rTor().data());

         // Set poloidal component to zero
         this->setZero(this->mpMagB->rOc().rPerturbation().rPol().data());

      } else if(TSimTraits::NeedMagnetic && TFilter::ReadMagnetic == StateFileDefs::PoloidalOnly)
      {
         // Read poloidal component
         this->readPolField(StateFileDefs::MAGNETICTAG, this->mpMagB->rOc().rPerturbation().rPol().data());

         // Set toroidal component to zero
         this->setZero(this->mpMagB->rOc().rPerturbation().rTor().data());
      }
      
      // Read requested velocity coefficients
      if(TSimTraits::NeedVelocity && TFilter::ReadVelocity == StateFileDefs::FullField)
      {
         this->readTorPolField(StateFileDefs::VELOCITYTAG, this->mpVelV->rOc().rPerturbation().rTor().data(), this->mpVelV->rOc().rPerturbation().rPol().data());
      } else if(TSimTraits::NeedVelocity && TFilter::ReadVelocity == StateFileDefs::ToroidalOnly)
      {
         // Read toroidal component
         this->readTorField(StateFileDefs::VELOCITYTAG, this->mpVelV->rOc().rPerturbation().rTor().data());

         // Set poloidal component to zero
         this->setZero(this->mpVelV->rOc().rPerturbation().rPol().data());
      } else if(TSimTraits::NeedVelocity && TFilter::ReadVelocity == StateFileDefs::PoloidalOnly)
      {
         // Read poloidal component
         this->readPolField(StateFileDefs::VELOCITYTAG, this->mpVelV->rOc().rPerturbation().rPol().data());

         // Set toroidal component to zero
         this->setZero(this->mpVelV->rOc().rPerturbation().rTor().data());
      }

   }

}

#endif // STATEFILEREADER_HPP
