/** \file TimestepControlBase.cpp
 *  \brief Source of the implementation of a simple timestep control
 */

// System includes
//
#include <limits>
#include <iostream>

// Configuration includes
//
#include "Config/Parallelisation.h"

// External includes
//

// Class include
//
#include "Timestepping/TimestepControlBase.hpp"

// Project includes
//
#include "Timestepping/TimestepConfig.hpp"

namespace EPMDynamo {

   TimestepControlBase::TimestepControlBase(TimestepParameters &tsParams, const EquationParameters &eqParams, TimestepCtrlTypes ctrlType, int order)
      : mKeepRunning(true), mNeedInit(true), mError(TimestepConfig::TIMESTEP_ERROR_EPSILON), mOldError(TimestepConfig::TIMESTEP_ERROR_EPSILON), mrTSParams(tsParams), mrEqParams(eqParams), mController(ctrlType, order, tsParams), mCtrlPI42(PI42Ctrl, order, tsParams), mCtrlH211B(H211BCtrl, order, tsParams)
   {
      // Set the two global CFL conditions
      this->rTSParams().setInertialCFL(this->eqParams().inertialCFL());
      this->rTSParams().setTorsionalCFL(this->eqParams().torsionalCFL());
   }

   void TimestepControlBase::resetError()
   {
      // Store error as old value
      this->mOldError = this->mError;
      
      // Reset error to zero
      this->mError = 0.0;
   }

   void TimestepControlBase::updateRTPCFLTimestep(const RTPField &velV)
   {
      if(this->rTSParams().isNextStep())
      {
         /////////////////////////////////////////////////////
         // ADDITIONAL ELEMENTS TEMPORARY REQUIRED
         EPMFloat radVPos;
         EPMFloat hozVPos;
         /////////////////////////////////////////////////////


         // Get truncation information
         int r0 = velV.trunc()->local()->rtp()->r0();
         int nR = velV.nR();
         int nRSim = velV.trunc()->sim()->rad()->nR();

         // Get radii
         Array radii = velV.trunc()->sim()->rad()->radGrid();
         // Get radial dependent sqrt(l(l+)) factor
         Array radSll = velV.trunc()->sim()->rad()->radSll();

         // Initialise to highest representable value
         EPMFloat minVRad = std::numeric_limits<EPMFloat>::max();
         EPMFloat minVHoz = std::numeric_limits<EPMFloat>::max();

         // Local CFL Condition
         EPMFloat dr;
         EPMFloat maxVel = 0.0;
         Matrix   p;
         int n_;
         for(int n=0; n < nR; ++n)
         {
            n_ = n + r0;
            if(n_ == 0)
            {
               dr = radii(1)-radii(0);
            }
            else if(n_ == nRSim - 1)
            {
               dr = 1.0 - radii(n_);
            } else
            {
               dr = std::min(radii(n_) - radii(n_-1), radii(n_+1) - radii(n_));
            }

            // Radial Velocity
            maxVel = velV.r().shell(n).array().abs().maxCoeff();

            // Update timestep
            if(maxVel != 0.0)
            {
               if(dr/maxVel < minVRad)
               {
                  minVRad = dr/maxVel;
                  radVPos = radii(n_);
               }
            }

            // Angular "Velocity"
            dr = radii(n_)/radSll(n_);
            maxVel = (velV.theta().shell(n).array().square() + velV.phi().shell(n).array().square()).array().sqrt().maxCoeff();

            // Update timestep
            if(maxVel != 0.0)
            {
               if(dr/maxVel < minVHoz)
               {
                  minVHoz = dr/maxVel;
                  hozVPos = radii(n_);
               }
            }
         }

         // Set the local CFL conditions from the velocity field
         this->rTSParams().setRadRTPVCFL(minVRad, radVPos);
         this->rTSParams().setHozRTPVCFL(minVHoz, hozVPos);
      }
   }

   void TimestepControlBase::updateRTPCFLTimestep(const RTPField &magB, const RTPField &velV)
   {
      if(this->rTSParams().isNextStep())
      {
         /////////////////////////////////////////////////////
         // ADDITIONAL ELEMENTS TEMPORARY REQUIRED
         EPMFloat radVBPos;
         EPMFloat hozVBPos;
         ////////////////////////////////////////////////////


         // Get truncation information
         int r0 = magB.trunc()->local()->rtp()->r0();
         int nR = magB.nR();
         int nRSim = magB.trunc()->sim()->rad()->nR();

         // Get radii
         Array radii = magB.trunc()->sim()->rad()->radGrid();
         // Get radial dependent sqrt(l(l+)) factor
         Array radSll = magB.trunc()->sim()->rad()->radSll();

         // Initialise to highest representable value
         EPMFloat minVBRad = std::numeric_limits<EPMFloat>::max();
         EPMFloat minVBHoz = std::numeric_limits<EPMFloat>::max();

         // Local CFL Condition
         EPMFloat d;
         EPMFloat dr;
         EPMFloat maxVel = 0.0;
         Matrix   p;
         int n_;
         for(int n=0; n < nR; ++n)
         {
            n_ = n + r0;
            if(n_ == 0)
            {
               dr = radii(1)-radii(0);
            }
            else if(n_ == nRSim - 1)
            {
               dr = 1.0 - radii(n_);
            } else
            {
               dr = std::min(radii(n_) - radii(n_-1), radii(n_+1) - radii(n_));
            }
            d = this->eqParams().alfvenDamping(dr);

            // Radial "Velocity"
            p = magB.r().shell(n).array().square()*this->eqParams().alfvenFactor();
            maxVel = (p.array()/(p.array() + d).array().sqrt() + velV.r().shell(n).array().abs()).maxCoeff();

            // Update timestep
            if(maxVel != 0.0)
            {
               if(dr/maxVel < minVBRad)
               {
                  minVBRad = dr/maxVel;
                  radVBPos = radii(n_);
               }
            }

            // Angular "Velocity"
            dr = radii(n_)/radSll(n_);
            d = this->eqParams().alfvenDamping(dr);
            p = (magB.theta().shell(n).array().square() + magB.phi().shell(n).array().square())*this->eqParams().alfvenFactor();
            maxVel = (p.array()/(p.array() + d).array().sqrt() + (velV.theta().shell(n).array().square() + velV.phi().shell(n).array().square()).array().sqrt()).maxCoeff();

            // Update timestep
            if(maxVel != 0.0)
            {
               if(dr/maxVel < minVBHoz)
               {
                  minVBHoz = dr/maxVel;
                  hozVBPos = radii(n_);
               }
            }
         }

         // Store the resulting CFL conditions
         this->rTSParams().setRadRTPVBCFL(minVBRad, radVBPos);
         this->rTSParams().setHozRTPVBCFL(minVBHoz, hozVBPos);

         /////////////////////////////////////////////////////
         // STORE TEMPORARY SOME ADDITIONAL INFORMATION
         this->updateRTPCFLTimestep(velV);
         ///////////////////////////////////////////////////
      }
   }

   void TimestepControlBase::updateSpecCFLTimestep(const SpectralSHScalar& velT, const SpectralSHScalar& velP)
   {
      if(this->rTSParams().isNextStep())
      {
         /////////////////////////////////////////////////////
         // STORE TEMPORARY SOME ADDITIONAL INFORMATION
         EPMFloat polRadVPos;
         int polRadVDeg;
         EPMFloat torHozVPos;
         int torHozVDeg;
         EPMFloat polHozVPos;
         int polHozVDeg;
         /////////////////////////////////////////////////////



         // Get truncation information
         int nR = velT.trunc()->sim()->rad()->nR();
         ArrayI ls = velT.trunc()->local()->spec()->lArray();
         int nL = ls.size();

         /////////////////////////////////////////////////////
         Array startIdx(nL);
         startIdx << 0,1,2,2,3,4,4,5,5,6,6,7,7,7,8,8,8,9,9,10,10,10,10,11,11,11,12,12,12,12,13,13;
         /////////////////////////////////////////////////////

         // Get radii
         Array radii = velT.trunc()->sim()->rad()->radGrid();

         // Compute the L2 norms by harmonic degree
         Array velTNorm = velT.l2NormByL();
         Array velPNorm = velP.l2NormByL();

         // Initialise to highest representable value
         EPMFloat minVPolRad = std::numeric_limits<EPMFloat>::max();
         EPMFloat minVTorHoz = std::numeric_limits<EPMFloat>::max();
         EPMFloat minVPolHoz = std::numeric_limits<EPMFloat>::max();

         // Local CFL Condition
         EPMFloat dr;
         int n_;
         int l;
         for(int i = 0; i < nL; i++)
         {
            // Get the harmonic degree l
            l = ls(i);

            // Local CFL Condition
            n_ = startIdx(l);

            for(;n_ < nR; n_++)
            {
               if(n_ == 0)
               {
                  dr = radii(1)-radii(0);
               }
               else if(n_ == nR - 1)
               {
                  dr = 1.0 - radii(n_);
               } else
               {
                  dr = std::min(radii(n_) - radii(n_-1), radii(n_+1) - radii(n_));
               }

               if(dr/velPNorm(l) < minVPolRad)
               {
                  minVPolRad = dr/velPNorm(l);
                  polRadVPos = radii(n_);
                  polRadVDeg = l;
               }

               dr = radii(n_)/std::sqrt(l*(l+1));

               if(dr/velTNorm(l) < minVTorHoz)
               {
                  minVTorHoz = dr/velTNorm(l);
                  torHozVPos = radii(n_);
                  torHozVDeg = l;
               }

               if(dr/velPNorm(l) < minVPolHoz)
               {
                  minVPolHoz = dr/velPNorm(l);
                  polHozVPos = radii(n_);
                  polHozVDeg = l;
               }
            }
         }

         // Store local spectral CFL conditions
         this->rTSParams().setPolRadSpecVCFL(minVPolRad, polRadVPos, polRadVDeg);
         this->rTSParams().setTorHozSpecVCFL(minVTorHoz, torHozVPos, torHozVDeg);
         this->rTSParams().setPolHozSpecVCFL(minVPolHoz, polHozVPos, polHozVDeg);
      }
   }

   void TimestepControlBase::updateSpecCFLTimestep(const SpectralSHScalar& magT, const SpectralSHScalar& magP, const SpectralSHScalar& velT, const SpectralSHScalar& velP)
   {
      if(this->rTSParams().isNextStep())
      {
         /////////////////////////////////////////////////////
         // STORE TEMPORARY SOME ADDITIONAL INFORMATION
         EPMFloat polRadVBPos;
         int polRadVBDeg;
         EPMFloat torHozVBPos;
         int torHozVBDeg;
         EPMFloat polHozVBPos;
         int polHozVBDeg;
         /////////////////////////////////////////////////////
         
         // Get truncation information
         int nR = velT.trunc()->sim()->rad()->nR();
         ArrayI ls = velT.trunc()->local()->spec()->lArray();
         int nL = ls.size();

         // Initialise to highest representable value
         EPMFloat minVBPolRad = std::numeric_limits<EPMFloat>::max();
         EPMFloat minVBTorHoz = std::numeric_limits<EPMFloat>::max();
         EPMFloat minVBPolHoz = std::numeric_limits<EPMFloat>::max();

         // Get radii
         Array radii = velT.trunc()->sim()->rad()->radGrid();
         
         Array startIdx(nL);
         startIdx << 0,1,2,2,3,4,4,5,5,6,6,7,7,7,8,8,8,9,9,10,10,10,10,11,11,11,12,12,12,12,13,13;

         Array magTNorm = magT.l2NormByL();
         Array magPNorm = magP.l2NormByL();
         Array velTNorm = velT.l2NormByL();
         Array velPNorm = velP.l2NormByL();

         // Local CFL Condition
         EPMFloat dr;
         EPMFloat d;
         EPMFloat maxVel = 0.0;
         EPMFloat p;
         int n_;
         int l;
         for(int i = 0; i < nL; i++)
         {
            // Get the right harmonic degree l
            l = ls(i);

            // Local CFL Condition
            n_ = startIdx(l);

            // Loop over the remaining radial grid points
            for(;n_ < nR; n_++)
            {
               if(n_ == 0)
               {
                  dr = radii(1)-radii(0);
               }
               else if(n_ == nR - 1)
               {
                  dr = 1.0 - radii(n_);
               } else
               {
                  dr = std::min(radii(n_) - radii(n_-1), radii(n_+1) - radii(n_));
               }
               d = this->eqParams().alfvenDamping(dr);

               p = magPNorm(l)*magPNorm(l)*this->eqParams().alfvenFactor();
               maxVel = p/std::sqrt(p + d) + velPNorm(l);

               if(dr/maxVel < minVBPolRad)
               {
                  minVBPolRad = dr/maxVel;
                  polRadVBPos = radii(n_);
                  polRadVBDeg = l;
               }

               dr = radii(n_)/std::sqrt(l*(l+1));
               d = this->eqParams().alfvenDamping(dr);
               p = magTNorm(l)*magTNorm(l)*this->eqParams().alfvenFactor();
               maxVel = p/std::sqrt(p + d) + velTNorm(l);

               if(dr/maxVel < minVBTorHoz)
               {
                  minVBTorHoz = dr/maxVel;
                  torHozVBPos = radii(n_);
                  torHozVBDeg = l;
               }

               p = magPNorm(l)*magPNorm(l)*this->eqParams().alfvenFactor();
               maxVel = p/std::sqrt(p + d) + velPNorm(l);

               if(dr/maxVel < minVBPolHoz)
               {
                  minVBPolHoz = dr/maxVel;
                  polHozVBPos = radii(n_);
                  polHozVBDeg = l;
               }
            }
         }

         // Store local CFL conditions
         this->rTSParams().setPolRadSpecVBCFL(minVBPolRad, polRadVBPos, polRadVBDeg);
         this->rTSParams().setTorHozSpecVBCFL(minVBTorHoz, torHozVBPos, torHozVBDeg);
         this->rTSParams().setPolHozSpecVBCFL(minVBPolHoz, polHozVBPos, polHozVBDeg);

         /////////////////////////////////////////////////////
         // STORE TEMPORARY SOME ADDITIONAL INFORMATION
         this->updateSpecCFLTimestep(velT, velP);
         /////////////////////////////////////////////////////
      }
   }

   EPMFloat TimestepControlBase::getSimulationCFLCondition()
   {
      EPMFloat cfl = this->rTSParams().getCFL();

      // Get the "global" local minimum for MPI code
      #ifdef EPMDYNAMO_MPI
         MPI_Allreduce(MPI_IN_PLACE, &cfl, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         // STORE TEMPORARY SOME ADDITIONAL INFORMATION
         // Get global RTP CFL minima
         MPI_Allreduce(MPI_IN_PLACE, this->rTSParams().rRTPCFLs().data(), this->rTSParams().rtpCFLs().size(), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

         // Get global spectral CFL minima
         Matrix   test(2,this->rTSParams().specCFLs().size());
         test.row(0) = this->rTSParams().specCFLs();
         test.row(1) = this->rTSParams().specCFLPos();
         MPI_Allreduce(MPI_IN_PLACE, test.data(), this->rTSParams().specCFLs().size(), MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD);
         this->rTSParams().rSpecCFLPos() = test.row(1);

         test.row(0) = this->rTSParams().specCFLs();
         test.row(1) = this->rTSParams().specCFLDeg();
         MPI_Allreduce(MPI_IN_PLACE, test.data(), this->rTSParams().specCFLs().size(), MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD);
         this->rTSParams().rSpecCFLDeg() = test.row(1);
         this->rTSParams().rSpecCFLs() = test.row(0);

         // Get global errorr CFL minima
         MPI_Allreduce(MPI_IN_PLACE, this->rTSParams().rErrCFLs().data(), this->rTSParams().errCFLs().size(), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      #endif // EPMDYNAMO_MPI

      return cfl;
   }

   void TimestepControlBase::testInitialisation(EPMFloat& rDt)
   {
      // Initialisation timestep
      if(this->mNeedInit)
      {
         // Set simple starting timestep to mininal timestep
         if(rDt == TimestepConfig::MIN_TIMESTEP)
         {
            rDt = TimestepConfig::MIN_TIMESTEP*TimestepConfig::TIMESTEP_MAX_RATIO;
         } else
         {
            rDt = std::max(rDt/1.0e3, TimestepConfig::MIN_TIMESTEP*TimestepConfig::TIMESTEP_MAX_RATIO);
         }

         // Set flag to forbid second entry
         this->mNeedInit = false;

         std::cout << "Use Initialisation Timestep!" << std::endl;
      }
   }

   void TimestepControlBase::testMaximumTimestep(EPMFloat& rDt)
   {
      // Force maximum timestep
      if(rDt > TimestepConfig::MAXIMUM_TIMESTEP)
      {
         // Set timestep to maximum allowed timestep
         rDt = TimestepConfig::MAXIMUM_TIMESTEP;
      }
   }

   void TimestepControlBase::testCFLCondition(EPMFloat& rDt)
   {
      // Get the simulation wide CFL condition
      EPMFloat cfl = this->getSimulationCFLCondition();

      // Get timestep according to CFL condition
      if(cfl > 0.0)
      {
         rDt = std::min(cfl, rDt);
      }
   }

   void TimestepControlBase::useErrorCtrlTimestep()
   {
      // Store the values for three different controllers
      this->rTSParams().setErrorCFL(this->mController.nextTimestep(this->mError, this->mOldError), 0);
      this->rTSParams().setErrorCFL(this->mCtrlPI42.nextTimestep(this->mError, this->mOldError), 1);
      this->rTSParams().setErrorCFL(this->mCtrlH211B.nextTimestep(this->mError, this->mOldError), 2);
   }

   void TimestepControlBase::useCourantTimestep(EPMFloat& rDt)
   {
      // Include Courant number into timestep value
      rDt = rDt*TimestepConfig::COURANT_NUMBER;
   }

   void TimestepControlBase::useWindowedTimestep(EPMFloat dt)
   {
      // Don't make adaptative scheme to reactive. Allow a range of timestep values
      if(dt < this->rTSParams().dt()*(1.0-TimestepConfig::TIMESTEP_MARGIN) ||  dt > this->rTSParams().dt()*(1.0 + TimestepConfig::TIMESTEP_MARGIN))
      {
         // Store new timestep if different from old step
         if(dt != this->rTSParams().dt())
         {
            this->rTSParams().updateTimestep(dt);

            std::cout << "Updating Timestep! Now using dt=" << this->rTSParams().dt() << " at time t=" << this->rTSParams().time()<< std::endl;
         }
      }
   }

   void TimestepControlBase::useTimestep(EPMFloat dt)
   {
//      dt = 5e-3;
//      std::cout << "TIMESTEP HAS BEEN LOCKED TO FIXED VALUE dt=" << dt << std::endl;

      // Store new timestep if different from old step
      if(dt != this->rTSParams().dt())
      {
         this->rTSParams().updateTimestep(dt);

         std::cout << "Updating Timestep! Now using dt=" << this->rTSParams().dt() << " at time t=" << this->rTSParams().time()<< std::endl;
      }
   }

}
