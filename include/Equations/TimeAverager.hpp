/** \file TimeAverager.hpp
 *  \brief Implementation of the time averaged field computation
 */

#ifndef TIMEAVERAGER_HPP
#define TIMEAVERAGER_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Timestepping/TimestepParameters.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the time averaged field computation
    *
    * \tparam TSimTraits Traits of the implementation
    */
   template <typename TSimTraits> class TimeAverager
   {
      public:
         typedef EPMSHARED_PTR<typename TSimTraits::CodType>  SmartCod;
         typedef EPMSHARED_PTR<typename TSimTraits::VelType>  SmartVel;
         typedef EPMSHARED_PTR<typename TSimTraits::MagType>  SmartMag;

         /// Typedef from Simulation trait to local transform type
         typedef SimulationConfig::TransformType    TransformType;

         /**
         * @brief Constructor with all three fields
         *
         * @param codC Codensity variable
         * @param magB Magnetic variable
         * @param velV Velocity variable
         * @param tsParams Timestep parameters
         */
         TimeAverager(const typename TSimTraits::CodType &codC, const typename TSimTraits::MagType &magB, const typename TSimTraits::VelType &velV, TransformType &transform, const TimestepParameters &tsParams);

         /**
         * @brief Constructor for thermal convection
         *
         * @param codC Codensity variable
         * @param velV Velocity variable
         * @param tsParams Timestep parameters
         */
         TimeAverager(const typename TSimTraits::CodType &codC, const typename TSimTraits::VelType &velV, TransformType &transform, const TimestepParameters &tsParams);

         /**
         * @brief Constructor for magneto-convection
         *
         * @param magB Magnetic variable
         * @param velV Velocity variable
         * @param tsParams Timestep parameters
         */
         TimeAverager(const typename TSimTraits::MagType &magB, const typename TSimTraits::VelType &velV, TransformType &transform, const TimestepParameters &tsParams);

         /**
         * @brief Constructor with only codensity
         *
         * @param codC Codensity variable
         * @param tsParams Timestep parameters
         */
         TimeAverager(const typename TSimTraits::CodType &codC, TransformType &transform, const TimestepParameters &tsParams);

         /**
         * @brief Constructor with only magnetic field
         *
         * @param magB Magnetic variable
         * @param tsParams Timestep parameters
         */
         TimeAverager(const typename TSimTraits::MagType &magB, TransformType &transform, const TimestepParameters &tsParams);

         /**
         * @brief Constructor for cases with all three fields
         *
         * @param velV Velocity variable
         * @param tsParams Timestep parameters
         */
         TimeAverager(const typename TSimTraits::VelType &velV, TransformType &transform, const TimestepParameters &tsParams);

         /**
         * @brief Destructor
         */
         virtual ~TimeAverager() {};

         /**
          * @brief Initialise the averager
          */
         void init();

         /**
          * @brief Include current timestep into average
          */
         void timestep();

         /**
          * @brief Get Time averaged codensity field
          */
         SmartCod pAvgCodC();

         /**
          * @brief Get Time averaged velocity field
          */
         SmartVel pAvgVelV();

         /**
          * @brief Get Time averaged magnetic field
          */
         SmartMag pAvgMagB();
         
      protected:
         /**
          * @brief Pointer to the condensity variable
          */
         const typename TSimTraits::CodType*  mpCodC;

         /**
          * @brief Pointer to the magnetic variable
          */
         const typename TSimTraits::MagType*  mpMagB;

         /**
          * @brief Pointer to the velocity variable
          */
         const typename TSimTraits::VelType*  mpVelV;

      private:
         /**
          * @brief Storage for the timestep to use in integration
          */
         EPMFloat mTimestep;

         /**
          * @brief Elapsed time since start of integration
          */
         EPMFloat mTime;

         /**
          * @brief Spectral truncatio information
          */
         SmartTruncation mpTrunc;

         /**
          * @brief Reference to a TimestepParameters object to get access to time and timestep
          */
         const TimestepParameters&  mrTSParams;

         /**
          * @brief Storage for the averaged codensity field
          */
         SmartCod  mpAvgCodC;

         /**
          * @brief Storage for the averaged velocity field
          */
         SmartVel  mpAvgVelV;

         /**
          * @brief Storage for the averaged magnetic field
          */
         SmartMag  mpAvgMagB;

         /**
          * @brief Storage for the old codensity field
          */
         SmartCod  mpOldCodC;

         /**
          * @brief Storage for the old velocity field
          */
         SmartVel  mpOldVelV;

         /**
          * @brief Storage for the old magnetic field
          */
         SmartMag  mpOldMagB;

   };

   template <typename TSimTraits> typename TimeAverager<TSimTraits>::SmartCod TimeAverager<TSimTraits>::pAvgCodC()
   {
      return this->mpAvgCodC;
   }

   template <typename TSimTraits> typename TimeAverager<TSimTraits>::SmartVel TimeAverager<TSimTraits>::pAvgVelV()
   {
      return this->mpAvgVelV;
   }

   template <typename TSimTraits> typename TimeAverager<TSimTraits>::SmartMag TimeAverager<TSimTraits>::pAvgMagB()
   {
      return this->mpAvgMagB;
   }

   template <typename TSimTraits> TimeAverager<TSimTraits>::TimeAverager(const typename TSimTraits::CodType &codC, const typename TSimTraits::MagType &magB, const typename TSimTraits::VelType &velV, typename TimeAverager<TSimTraits>::TransformType &transform, const TimestepParameters &tsParams)
      : mTimestep(0.0), mTime(0.0), mpTrunc(codC.oc().trunc()), mpCodC(&codC), mpMagB(&magB), mpVelV(&velV), mrTSParams(tsParams)
   {
      // Initialise the codensity storage and average
      this->mpAvgCodC = SmartCod(new typename TSimTraits::CodType(this->mpTrunc, transform));
      this->mpAvgCodC->initialiseZeros();

      this->mpOldCodC = SmartCod(new typename TSimTraits::CodType(this->mpTrunc, transform));
      this->mpOldCodC->initialiseZeros();

      // Initialise the velocity storage and average
      this->mpAvgVelV = SmartVel(new typename TSimTraits::VelType(this->mpTrunc, transform));
      this->mpAvgVelV->initialiseZeros();

      this->mpOldVelV = SmartVel(new typename TSimTraits::VelType(this->mpTrunc, transform));
      this->mpOldVelV->initialiseZeros();

      // Initialise the magnetic storage and average
      this->mpAvgMagB = SmartMag(new typename TSimTraits::MagType(this->mpTrunc, transform));
      this->mpAvgMagB->initialiseZeros();

      this->mpOldMagB = SmartMag(new typename TSimTraits::MagType(this->mpTrunc, transform));
      this->mpOldMagB->initialiseZeros();
   }

   template <typename TSimTraits> TimeAverager<TSimTraits>::TimeAverager(const typename TSimTraits::CodType &codC, const typename TSimTraits::VelType &velV, typename TimeAverager<TSimTraits>::TransformType &transform, const TimestepParameters &tsParams)
      : mTimestep(0.0), mTime(0.0), mpTrunc(codC.oc().trunc()), mpCodC(&codC), mpMagB(NULL), mpVelV(&velV), mrTSParams(tsParams)
   {
      // Initialise the codensity storage and average
      this->mpAvgCodC = SmartCod(new typename TSimTraits::CodType(this->mpTrunc, transform));
      this->mpAvgCodC->initialiseZeros();

      this->mpOldCodC = SmartCod(new typename TSimTraits::CodType(this->mpTrunc, transform));
      this->mpOldCodC->initialiseZeros();

      // Initialise the velocity storage and average
      this->mpAvgVelV = SmartVel(new typename TSimTraits::VelType(this->mpTrunc, transform));
      this->mpAvgVelV->initialiseZeros();

      this->mpOldVelV = SmartVel(new typename TSimTraits::VelType(this->mpTrunc, transform));
      this->mpOldVelV->initialiseZeros();
   }

   template <typename TSimTraits> TimeAverager<TSimTraits>::TimeAverager(const typename TSimTraits::MagType &magB, const typename TSimTraits::VelType &velV, typename TimeAverager<TSimTraits>::TransformType &transform, const TimestepParameters &tsParams)
      : mTimestep(0.0), mTime(0.0), mpTrunc(magB.oc().trunc()), mpCodC(NULL), mpMagB(&magB), mpVelV(&velV), mrTSParams(tsParams)
   {
      // Initialise the velocity storage and average
      this->mpAvgVelV = SmartVel(new typename TSimTraits::VelType(this->mpTrunc, transform));
      this->mpAvgVelV->initialiseZeros();

      this->mpOldVelV = SmartVel(new typename TSimTraits::VelType(this->mpTrunc, transform));
      this->mpOldVelV->initialiseZeros();

      // Initialise the magnetic storage and average
      this->mpAvgMagB = SmartMag(new typename TSimTraits::MagType(this->mpTrunc, transform));
      this->mpAvgMagB->initialiseZeros();

      this->mpOldMagB = SmartMag(new typename TSimTraits::MagType(this->mpTrunc, transform));
      this->mpOldMagB->initialiseZeros();
   }

   template <typename TSimTraits> TimeAverager<TSimTraits>::TimeAverager(const typename TSimTraits::CodType &codC, typename TimeAverager<TSimTraits>::TransformType &transform, const TimestepParameters &tsParams)
      : mTimestep(0.0), mTime(0.0), mpTrunc(codC.oc().trunc()), mpCodC(&codC), mpMagB(NULL), mpVelV(NULL),  mrTSParams(tsParams)
   {
      // Initialise the codensity storage and average
      this->mpAvgCodC = SmartCod(new typename TSimTraits::CodType(this->mpTrunc, transform));
      this->mpAvgCodC->initialiseZeros();

      this->mpOldCodC = SmartCod(new typename TSimTraits::CodType(this->mpTrunc, transform));
      this->mpOldCodC->initialiseZeros();
   }

   template <typename TSimTraits> TimeAverager<TSimTraits>::TimeAverager(const typename TSimTraits::MagType &magB, typename TimeAverager<TSimTraits>::TransformType &transform, const TimestepParameters &tsParams)
      : mTimestep(0.0), mTime(0.0), mpTrunc(magB.oc().trunc()), mpCodC(NULL), mpMagB(&magB), mpVelV(NULL), mrTSParams(tsParams)
   {
      // Initialise the magnetic storage and average
      this->mpAvgMagB = SmartMag(new typename TSimTraits::MagType(this->mpTrunc, transform));
      this->mpAvgMagB->initialiseZeros();

      this->mpOldMagB = SmartMag(new typename TSimTraits::MagType(this->mpTrunc, transform));
      this->mpOldMagB->initialiseZeros();
   }

   template <typename TSimTraits> TimeAverager<TSimTraits>::TimeAverager(const typename TSimTraits::VelType &velV, typename TimeAverager<TSimTraits>::TransformType &transform, const TimestepParameters &tsParams)
      : mTimestep(0.0), mTime(0.0), mpTrunc(velV.oc().trunc()), mpCodC(NULL), mpMagB(NULL), mpVelV(&velV), mrTSParams(tsParams)
   {
      // Initialise the velocity storage and average
      this->mpAvgVelV = SmartVel(new typename TSimTraits::VelType(this->mpTrunc, transform));
      this->mpAvgVelV->initialiseZeros();

      this->mpOldVelV = SmartVel(new typename TSimTraits::VelType(this->mpTrunc, transform));
      this->mpOldVelV->initialiseZeros();
   }

   template <typename TSimTraits> void TimeAverager<TSimTraits>::init()
   {
   }

   template <typename TSimTraits> void TimeAverager<TSimTraits>::timestep()
   {
      // Make sure we are not in some intermediate stage
      if(this->mrTSParams.isNextStep())
      {
         // Get number of harmonic degrees
         int nL;
         int l0;

         EPMFloat factor = this->mTime/(this->mTime+this->mTimestep);
         EPMFloat weight = 0.5 * this->mTimestep / (this->mTime + this->mTimestep);

         if(this->mTime + this->mTimestep == 0.0)
         {
            factor = 0.0;
            weight = 0.0;
         }

         // Update average of the codensity storage and average
         if(this->mpCodC != NULL)
         {
            nL = this->mpCodC->oc().perturbation().trunc()->local()->spec()->nL();
            l0 = this->mpCodC->oc().perturbation().minL();

            for(int l = l0; l < nL; l++)
            {
               this->mpAvgCodC->rOc().rPerturbation().rLShell(l) *= factor;

               this->mpAvgCodC->rOc().rPerturbation().rLShell(l) += (this->mpOldCodC->oc().perturbation().lshell(l) + this->mpCodC->oc().perturbation().lshell(l))*weight;

               this->mpOldCodC->rOc().rPerturbation().rLShell(l) = this->mpCodC->oc().perturbation().lshell(l);
            }
         }

         // Update average of the velocity storage and average
         if(this->mpVelV != NULL)
         {
            nL = this->mpVelV->oc().perturbation().tor().trunc()->local()->spec()->nL();
            l0 = this->mpVelV->oc().perturbation().tor().minL();

            for(int l = l0; l < nL; l++)
            {
               this->mpAvgVelV->rOc().rPerturbation().rTor().rLShell(l) *= factor;

               this->mpAvgVelV->rOc().rPerturbation().rTor().rLShell(l) += (this->mpOldVelV->oc().perturbation().tor().lshell(l) + this->mpVelV->oc().perturbation().tor().lshell(l))*weight;

               this->mpAvgVelV->rOc().rPerturbation().rPol().rLShell(l) *= factor;

               this->mpAvgVelV->rOc().rPerturbation().rPol().rLShell(l) += (this->mpOldVelV->oc().perturbation().pol().lshell(l) + this->mpVelV->oc().perturbation().pol().lshell(l))*weight;

               this->mpOldVelV->rOc().rPerturbation().rTor().rLShell(l) = this->mpVelV->oc().perturbation().tor().lshell(l);

               this->mpOldVelV->rOc().rPerturbation().rPol().rLShell(l) = this->mpVelV->oc().perturbation().pol().lshell(l);
            }
         }

         // Update average of the magnetic storage and average
         if(this->mpMagB != NULL)
         {
            nL = this->mpMagB->oc().perturbation().tor().trunc()->local()->spec()->nL();
            l0 = this->mpMagB->oc().perturbation().tor().minL();

            for(int l = l0; l < nL; l++)
            {
               this->mpAvgMagB->rOc().rPerturbation().rTor().rLShell(l) *= factor;

               this->mpAvgMagB->rOc().rPerturbation().rTor().rLShell(l) += (this->mpOldMagB->oc().perturbation().tor().lshell(l) + this->mpMagB->oc().perturbation().tor().lshell(l))*weight;

               this->mpAvgMagB->rOc().rPerturbation().rPol().rLShell(l) *= factor;

               this->mpAvgMagB->rOc().rPerturbation().rPol().rLShell(l) += (this->mpOldMagB->oc().perturbation().pol().lshell(l) + this->mpMagB->oc().perturbation().pol().lshell(l))*weight;

               this->mpOldMagB->rOc().rPerturbation().rTor().rLShell(l) = this->mpMagB->oc().perturbation().tor().lshell(l);

               this->mpOldMagB->rOc().rPerturbation().rPol().rLShell(l) = this->mpMagB->oc().perturbation().pol().lshell(l);
            }
         }

         this->mTime += this->mTimestep;

         this->mTimestep = this->mrTSParams.dt();
      }
   }

}

#endif // TIMEAVERAGER_HPP
