/** \file TorPolSHTransform.hpp
 *  \brief Implementation of the Toroidal/Poloidal spherical harmonics transforms
 */

#ifndef TORPOLSHTRANSFORME_HPP
#define TORPOLSHTRANSFORME_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "GeneralScalars/RTPScalar.hpp"
#include "GeneralFields/RTPField.hpp"
#include "Polynomials/AssocLegendreBasis.hpp"
#include "Polynomials/PolyOps/SetProduct.hpp"
#include "Polynomials/PolyOps/AddProduct.hpp"

namespace EPMDynamo {

   /**
    * @brief Implements the real space to FDSH representation transforms
    *
    * This transforms are required whenever the code is going to be an FDSH code or 
    * a full spectral code. This methods allow to perform the transformations between
    * real space and the Spherical harmonics expansion space. The radial behavior is not 
    * explicitly touched. (it is influenced in the sence that the spherical harmonics
    * expansion puts some contrains on the radial behavior for each harmonic mode.)
    *
    * \tparam TSHTraits Traits for the spherical harmonics transforms
    */
   template <typename TSHTraits> class TorPolSHTransform
   {
      public:
         /// Typedef for the Spherical harmonics scalar data type
         typedef typename TSHTraits::SHScalarType    SHScalarType;

         /// Typedef for the FFT forward data type
         typedef typename TSHTraits::FFTForwardType  FFTForwardType;

         /// Typedef for the FFT backward data type
         typedef typename TSHTraits::FFTBackwardType FFTBackwardType;

         /// Typedef for the DataManipulator data type
         typedef typename TSHTraits::DataManipulatorType   DataManipulatorType;

         /// Typedef for the Legendre Transform data type
         typedef typename TSHTraits::LegendreTransformType LegendreTransformType;

         /// Typedef for the FFT Transform data type
         typedef typename TSHTraits::FFTransformType   FFTransformType;

         /**
         * @brief Constructs the SH transform object
         *
         * This constructor will build up the two transforms that are required do do 
         * a spherical harmonics expansion: an FFT and a Legendre transform.
         *
         * @param pTrunc   Truncation information
         */
         TorPolSHTransform(const SmartTruncation pTrunc);

         /**
         * @brief Simple empty destructor
         */
         virtual ~TorPolSHTransform() {};

         /**
          * @brief Transform RTP scalar into its SH decomposition
          *
          * @param rSHValues  Output spherical harmonics coefficients
          * @param rtpValues  Input real space values
          */
         void transformRTP2Spec(SHScalarType   &rSHValues, const RTPScalar   &rtpValues);

         /**
          * @brief Transform RTP scalar into its SH decomposition
          *
          * This routine does only a partial transform used in a two step transform
          *
          * @param rtpValues  Input real space values
          */
         void transformRTP2Spec(const RTPScalar   &rtpValues);

         /**
          * @brief Transform an SH decomposition into RTP scalar
          *
          * @param rRTPValues Output real space values
          * @param shValues   Input spherical harmonics coefficients
          */
         void transformSpec2RTP(RTPScalar   &rRTPValues, const SHScalarType   &shValues);

         /**
          * @brief Transform an SH decomposition into RTP scalar
          *
          * This routine doesn only a partial transform used in a two step transform
          *
          * @param rRTPValues Output real space values
          */
         void transformSpec2RTP(RTPScalar   &rRTPValues);

         /**
          * @brief Get basis of the horizontal transforms
          */
         AssocLegendreBasis&  hozBasis();

         /**
          * @brief Get basis of the horizontal transforms
          */
         const AssocLegendreBasis&  hozBasis() const;

         /**
          * @brief Get the data manipulator
          */
         DataManipulatorType&  shManipulator();

      protected:

         /**
          * @brief Transform (R)TP to QST decomposition: Theta part of ST components
          *
          * Both the S component and the T component depend on the Theta component.
          *
          * @param rSComp  S component of QST decomposition
          * @param rTComp  T component of QST decomposition
          * @param thetaValues  Input real space theta values
          */
         void transformTP2STPartT(SHScalarType   &rSComp, SHScalarType   &rTComp, const RTPScalar   &thetaValues);

         /**
          * @brief Transform (R)TP to QST decomposition: Theta part of ST components
          *
          * This routine doesn only a partial transform used in a two step transform
          * Both the S component and the T component depend on the Theta component.
          *
          * @param thetaValues  Input real space theta values
          */
         void transformTP2STPartT(const RTPScalar   &thetaValues);

         /**
          * @brief Transform (R)TP to QST decomposition: Phi part of ST components
          *
          * This routine doesn only a partial transform used in a two step transform
          * Both the S component and the T component depend on the Phi component.
          *
          * @param rSComp  S component of QST decomposition
          * @param rTComp  T component of QST decomposition
          * @param phiValues  Input real space phi values
          */
         void transformTP2STPartP(SHScalarType   &rSComp, SHScalarType   &rTComp, const RTPScalar  &phiValues);

         /**
          * @brief Transform (R)TP to QST decomposition: Phi part of ST components
          *
          * This routine doesn only a partial transform used in a two step transform
          * Both the S component and the T component depend on the Phi component.
          *
          * @param phiValues  Input real space phi values
          */
         void transformTP2STPartP(const RTPScalar   &phiValues);

         /**
          * @brief Transform (Q)ST decomposition to RTP: Get Theta Component
          *
          * The Theta component is obtaind by adding the contribution from the S and Q
          * components
          *
          * @param rThetaValues Output real space theta values
          * @param sComp  Input S component of QST field
          * @param tComp  Input T component of QST field
          */
         void transformST2Theta(RTPScalar  &rThetaValues, const SHScalarType &sComp, const SHScalarType &tComp);

         /**
          * @brief Transform (Q)ST decomposition to RTP: Get Theta Component
          *
          * This routine doesn only a partial transform used in a two step transform
          * The Theta component is obtaind by adding the contribution from the S and Q
          * components
          *
          * @param rThetaValues Output real space theta values
          */
         void transformST2Theta(RTPScalar  &rThetaValues);
         
         /**
          * @brief Transform (Q)ST decomposition to RTP: Get Phi Component
          *
          * The Theta component is obtaind by adding the contribution from the S and Q
          * components
          *
          * @param rPhiValues Output real space theta values
          * @param sComp  Input S component of QST field
          * @param tComp  Input T component of QST field
          */
         void transformST2Phi(RTPScalar  &rPhiValues, const SHScalarType &sComp, const SHScalarType &tComp);
         
         /**
          * @brief Transform (Q)ST decomposition to RTP: Get Phi Component
          *
          * This routine doesn only a partial transform used in a two step transform
          * The Theta component is obtaind by adding the contribution from the S and Q
          * components
          *
          * @param rPhiValues Output real space theta values
          */
         void transformST2Phi(RTPScalar  &rPhiValues);

         /**
          * @brief Transform TP Components into T component
          *
          * @param rTComp T Component from QST decomposition
          * @param thetaComp  R Component from RTP decomposition
          * @param phiComp  R Component from RTP decomposition
          */
         void transformTP2T(SHScalarType &rTComp, const RTPScalar &thetaComp, const RTPScalar &phiComp);

         /**
          * @brief Transform TP Components into T component
          *
          * This routine doesn only a partial transform used in a two step transform
          *
          * @param thetaComp  R Component from RTP decomposition
          * @param phiComp  R Component from RTP decomposition
          */
         void transformTP2T(const RTPScalar &thetaComp, const RTPScalar &phiComp);

         /**
          * @brief Transform Scalar to Theta and Phi component of scalar
          *
          * @param rThComp Theta Component from QST gradient
          * @param rPhComp  Phi Component from RTP gradient
          * @param scalar  Spectral scalar already including r dependy computation
          */
         void transformSpec2GradTP(RTPScalar &rThComp, RTPScalar &rPhComp, const SHScalarType &scalar);

         /**
          * @brief Transform Scalar to Theta and Phi component of scalar
          *
          * This routine doesn only a partial transform used in a two step transform
          *
          * @param rThComp Theta Component from QST gradient
          * @param rPhComp  Phi Component from RTP gradient
          */
         void transformSpec2GradTP(RTPScalar &rThComp, RTPScalar &rPhComp);

      private:

         /**
          * @brief Is the m=0 modes stored on this cpu?
          */
         bool  mHasM0;

         /**
          * @brief FFTransform object for the FFT part
          */
         FFTransformType   mFFTransform;

         /**
          * @brief LegendreTransform object for the Legendre part
          */
         LegendreTransformType   mLegTransform;

         /**
          * @brief Data manipulator
          */
         DataManipulatorType   mSHManipulator;

         /**
          * @brief Does CPU hold m=0 mode?
          */
         bool hasM0() const;

         /**
          * @brief Get the FFT transform
          */
         FFTransformType&  fft();

         /**
          * @brief Get the legendre transform
          */
         LegendreTransformType&  legT();

         /**
          * @brief Send part of transform of RTP scalar into its SH decomposition
          *
          * @param rtpValues  Input real space values
          */
         void send_RTP2Spec(const RTPScalar   &rtpValues);

         /**
          * @brief Receive part of the transform of a SH decomposition into RTP scalar
          *
          * @param rRTPValues Output real space values
          */
         void recv_Spec2RTP(RTPScalar   &rRTPValues);

         /**
          * @brief Receive part of transform Scalar to Theta and Phi component of scalar
          *
          * @param rThComp Theta Component from QST gradient
          * @param rPhComp  Phi Component from RTP gradient
          */
         void recv_Spec2GradTP(RTPScalar &rThComp, RTPScalar &rPhComp);

         /**
          * @brief Send part of transform TP Components into T component
          *
          * @param thetaComp  R Component from RTP decomposition
          * @param phiComp  R Component from RTP decomposition
          */
         void send_TP2T(const RTPScalar &thetaComp, const RTPScalar &phiComp);

         /**
          * @brief Receive part of transform (Q)ST decomposition to RTP: Get Phi Component
          *
          * The Theta component is obtaind by adding the contribution from the S and Q
          * components
          *
          * @param rPhiValues Output real space theta values
          */
         void recv_ST2Phi(RTPScalar  &rPhiValues);

         /**
          * @brief Receive part of transform (Q)ST decomposition to RTP: Get Theta Component
          *
          * The Theta component is obtaind by adding the contribution from the S and Q
          * components
          *
          * @param rThetaValues Output real space theta values
          */
         void recv_ST2Theta(RTPScalar  &rThetaValues);

         /**
          * @brief Send part of transform (R)TP to QST decomposition: Phi part of ST components
          *
          * Both the S component and the T component depend on the Phi component.
          *
          * @param phiValues  Input real space phi values
          */
         void send_TP2STPartP(const RTPScalar   &phiValues);

         /**
          * @brief Send part of transform (R)TP to QST decomposition: Theta part of ST components
          *
          * Both the S component and the T component depend on the Theta component.
          *
          * @param thetaValues  Input real space theta values
          */
         void send_TP2STPartT(const RTPScalar   &thetaValues);
   };

   template <typename TSHTraits> inline bool TorPolSHTransform<TSHTraits>::hasM0() const
   {
      return this->mHasM0;
   }

   template <typename TSHTraits> inline AssocLegendreBasis&  TorPolSHTransform<TSHTraits>::hozBasis()
   {
      return this->legT().basis();
   }

   template <typename TSHTraits> inline const AssocLegendreBasis&  TorPolSHTransform<TSHTraits>::hozBasis() const
   {
      return this->mLegTransform.basis();
   }

   template <typename TSHTraits> inline typename TorPolSHTransform<TSHTraits>::DataManipulatorType& TorPolSHTransform<TSHTraits>::shManipulator()
   {
      return this->mSHManipulator;
   }

   template <typename TSHTraits> inline typename TorPolSHTransform<TSHTraits>::FFTransformType& TorPolSHTransform<TSHTraits>::fft()
   {
      return this->mFFTransform;
   }

   template <typename TSHTraits> inline typename TorPolSHTransform<TSHTraits>::LegendreTransformType& TorPolSHTransform<TSHTraits>::legT()
   {
      return this->mLegTransform;
   }

   template <typename TSHTraits> TorPolSHTransform<TSHTraits>::TorPolSHTransform(const SmartTruncation pTrunc)
      : mHasM0(false), mFFTransform(pTrunc), mLegTransform(pTrunc), mSHManipulator(pTrunc)
   {
      if(pTrunc->local()->fdsh()->mArray()(0) == 0)
      {
         this->mHasM0 = true;
      }
   }

   template <typename TSHTraits> void TorPolSHTransform<TSHTraits>::transformRTP2Spec(const RTPScalar &rtpValues)
   {
      this->shManipulator().addPacks(1);

      this->send_RTP2Spec(rtpValues);
   }

   template <typename TSHTraits> void TorPolSHTransform<TSHTraits>::send_RTP2Spec(const RTPScalar &rtpValues)
   {
      // Get temporary storage
      FFTForwardType  &rFTmp = this->shManipulator().provideFTmp();

      // Do the FFT transform of the real space values
      this->fft().forward(rFTmp, rtpValues);

      // Send data and release temporary storage
      this->shManipulator().send(rFTmp);
      this->shManipulator().releaseTmp(rFTmp);
   }

   template <typename TSHTraits> void TorPolSHTransform<TSHTraits>::transformRTP2Spec(typename TorPolSHTransform<TSHTraits>::SHScalarType &rSHValues, const RTPScalar &rtpValues)
   {
      if(this->shManipulator().isFirstEntry())
      {
         // add Packet size of transform
         this->shManipulator().addPacks(1);
      }

      if(this->shManipulator().atSendStage())
      {
         this->send_RTP2Spec(rtpValues);
      }

      if(this->shManipulator().atRecvStage())
      {
         // Get temporary storage and receive data
         FFTBackwardType &rMTmp = this->shManipulator().provideBTmp();
         this->shManipulator().receive(rMTmp);

         // Perform integration over Legendre polynomials
         this->legT().template multM<SetProduct>(rSHValues, rMTmp, &AssociatedLegendreOperator::intg);

         // Free temporary data
         this->shManipulator().freeTmp(rMTmp);

         // Set imaginary part of m=0 mode to zero
         //    Need to check if m=0 ist stored on this cpu
         if(this->hasM0())
         {
            int nR = rSHValues.nR(0);
            int nL = rSHValues.nL(0);
            for(int j=0; j < nR; ++j)
            {
               for(int i=0; i < nL; ++i)
               {
                  rSHValues.rMShell(0)(i,j).imag() = 0.0;
               }
            }
         }
      }
   }

   template <typename TSHTraits> void TorPolSHTransform<TSHTraits>::transformSpec2RTP(RTPScalar &rRTPValues)
   {
      this->recv_Spec2RTP(rRTPValues);
   }

   template <typename TSHTraits> void TorPolSHTransform<TSHTraits>::recv_Spec2RTP(RTPScalar &rRTPValues)
   {
      // Get temporary storage and receive data
      FFTForwardType &rFTmp = this->shManipulator().provideFTmp();
      rFTmp.doZeroPadding();
      this->shManipulator().receive(rFTmp);

      // Do the FFT of the complex space values
      this->fft().backward(rRTPValues, rFTmp);

      // Free temporary data
      this->shManipulator().freeTmp(rFTmp);
   }

   template <typename TSHTraits> void TorPolSHTransform<TSHTraits>::transformSpec2RTP(RTPScalar &rRTPValues, const typename TorPolSHTransform<TSHTraits>::SHScalarType &shValues)
   {
      if(this->shManipulator().isFirstEntry())
      {
         // add Packet size of transform
         this->shManipulator().addPacks(1);
      }

      if(this->shManipulator().atSendStage())
      {
         // add Packet size of transform
         FFTBackwardType &rMTmp = this->shManipulator().provideBTmp();

         // Perform summation over Legendre polynomials
         this->legT().template multM<SetProduct>(rMTmp, shValues, &AssociatedLegendreOperator::proj);

         // Send data and release temporary storage
         this->shManipulator().send(rMTmp);
         this->shManipulator().releaseTmp(rMTmp);
      }

      if(this->shManipulator().atRecvStage())
      {
         this->recv_Spec2RTP(rRTPValues);
      }
   }

   template <typename TSHTraits> void TorPolSHTransform<TSHTraits>::transformTP2STPartT(const RTPScalar &thetaValues)
   {   
         this->shManipulator().addPacks(1);

         this->send_TP2STPartT(thetaValues);
   }

   template <typename TSHTraits> void TorPolSHTransform<TSHTraits>::send_TP2STPartT(const RTPScalar &thetaValues)
   {   
         // Get temporary storage
         FFTForwardType  &rFTmp = this->shManipulator().provideFTmp();

         // Do the FFT transform of the real space values
         this->fft().forward(rFTmp, thetaValues);

         // Send data and release temporary storage
         this->shManipulator().send(rFTmp);
         this->shManipulator().releaseTmp(rFTmp);
   }

   template <typename TSHTraits> void TorPolSHTransform<TSHTraits>::transformTP2STPartT(typename TorPolSHTransform<TSHTraits>::SHScalarType &rSComp, typename TorPolSHTransform<TSHTraits>::SHScalarType &rTComp, const RTPScalar &thetaValues)
   {   
      if(this->shManipulator().isFirstEntry())
      {
         // add Packet size of transform
         this->shManipulator().addPacks(1);
      }

      if(this->shManipulator().atSendStage())
      {
         this->send_TP2STPartT(thetaValues);
      }

      if(this->shManipulator().atRecvStage())
      {
         // Get temporary storage and receive data
         FFTBackwardType &rMTmp = this->shManipulator().provideBTmp();
         this->shManipulator().receive(rMTmp);

         // Perform integration over Legendre polynomials: compute S component
         this->legT().template multM<SetProduct>(rSComp, rMTmp, &AssociatedLegendreOperator::intgTh2S);

         // Perform integration over Legendre polynomials: compute T component
         this->legT().template multM<SetProduct>(rTComp, rMTmp, &AssociatedLegendreOperator::intgTh2T);

         // Free temporary data
         this->shManipulator().freeTmp(rMTmp);

         // Set imaginary part of m=0 mode to zero
         //    Need to check if m=0 ist stored on this cpu
         if(this->hasM0())
         {
            int nR = rSComp.nR(0);
            int nL = rSComp.nL(0);
            for(int j=0; j < nR; ++j)
            {
               for(int i=0; i < nL; ++i)
               {
                  rSComp.rMShell(0)(i,j).imag() = 0.0;
                  rTComp.rMShell(0)(i,j).imag() = 0.0;
               }
            }
         }
      }
   }

   template <typename TSHTraits> void TorPolSHTransform<TSHTraits>::transformTP2STPartP(const RTPScalar &phiValues)
   {
      this->shManipulator().addPacks(1);

      this->send_TP2STPartP(phiValues);
   }

   template <typename TSHTraits> void TorPolSHTransform<TSHTraits>::send_TP2STPartP(const RTPScalar &phiValues)
   {
      // Get temporary storage
      FFTForwardType  &rFTmp = this->shManipulator().provideFTmp();

      // Do the FFT transform of the real space values
      this->fft().forward(rFTmp, phiValues);

      // Send data and release temporary storage
      this->shManipulator().send(rFTmp);
      this->shManipulator().releaseTmp(rFTmp);
   }

   template <typename TSHTraits> void TorPolSHTransform<TSHTraits>::transformTP2STPartP(typename TorPolSHTransform<TSHTraits>::SHScalarType &rSComp, typename TorPolSHTransform<TSHTraits>::SHScalarType &rTComp, const RTPScalar &phiValues)
   {
      if(this->shManipulator().isFirstEntry())
      {
         // add Packet size of transform
         this->shManipulator().addPacks(1);
      }

      if(this->shManipulator().atSendStage())
      {
         this->send_TP2STPartP(phiValues);
      }

      if(this->shManipulator().atRecvStage())
      {
         // Get temporary storage and receive data
         FFTBackwardType &rMTmp = this->shManipulator().provideBTmp();
         this->shManipulator().receive(rMTmp);

         // Perform integration over Legendre polynomials: compute S component
         this->legT().template multM<AddProduct>(rSComp, rMTmp, &AssociatedLegendreOperator::intgPh2S);

         // Perform integration over Legendre polynomials: compute T component
         this->legT().template multM<AddProduct>(rTComp, rMTmp, &AssociatedLegendreOperator::intgPh2T);

         // Free temporary data
         this->shManipulator().freeTmp(rMTmp);

         // Set imaginary part of m=0 mode to zero
         //    Need to check if m=0 ist stored on this cpu
         if(this->hasM0())
         {
            int nR = rSComp.nR(0);
            int nL = rSComp.nL(0);
            for(int j=0; j < nR; ++j)
            {
               for(int i=0; i < nL; ++i)
               {
                  rSComp.rMShell(0)(i,j).imag() = 0.0;
                  rTComp.rMShell(0)(i,j).imag() = 0.0;
               }
            }
         }
      }
   }

   template <typename TSHTraits> void TorPolSHTransform<TSHTraits>::transformST2Theta(RTPScalar &rThetaValues)
   {
      this->recv_ST2Theta(rThetaValues);
   }

   template <typename TSHTraits> void TorPolSHTransform<TSHTraits>::recv_ST2Theta(RTPScalar &rThetaValues)
   {
      // Get temporary storage and receive data
      FFTForwardType &rFTmp = this->shManipulator().provideFTmp();
      rFTmp.doZeroPadding();
      this->shManipulator().receive(rFTmp);

      // Do the FFT of the complex space values
      this->fft().backward(rThetaValues, rFTmp);

      // Free temporary data
      this->shManipulator().freeTmp(rFTmp);
   }

   template <typename TSHTraits> void TorPolSHTransform<TSHTraits>::transformST2Theta(RTPScalar &rThetaValues, const typename TorPolSHTransform<TSHTraits>::SHScalarType &sComp, const typename TorPolSHTransform<TSHTraits>::SHScalarType &tComp)
   {
      if(this->shManipulator().isFirstEntry())
      {
         // add Packet size of transform
         this->shManipulator().addPacks(1);
      }

      if(this->shManipulator().atSendStage())
      {
         // add Packet size of transform
         FFTBackwardType &rMTmp = this->shManipulator().provideBTmp();

         // Perform summation over Legendre polynomials: compute S contribution
         this->legT().template multM<SetProduct>(rMTmp, sComp, &AssociatedLegendreOperator::projS2Th);
         // Perform summation over Legendre polynomials: compute T contribution
         this->legT().template multM<AddProduct>(rMTmp, tComp, &AssociatedLegendreOperator::projT2Th);

         // Send data and release temporary storage
         this->shManipulator().send(rMTmp);
         this->shManipulator().releaseTmp(rMTmp);
      }

      if(this->shManipulator().atRecvStage())
      {
         this->recv_ST2Theta(rThetaValues);
      }
   }

   template <typename TSHTraits> void TorPolSHTransform<TSHTraits>::transformST2Phi(RTPScalar &rPhiValues)
   {
      this->recv_ST2Phi(rPhiValues);
   }

   template <typename TSHTraits> void TorPolSHTransform<TSHTraits>::recv_ST2Phi(RTPScalar &rPhiValues)
   {
      // Get temporary storage and receive data
      FFTForwardType &rFTmp = this->shManipulator().provideFTmp();
      rFTmp.doZeroPadding();
      this->shManipulator().receive(rFTmp);

      // Do the FFT of the complex space values
      this->fft().backward(rPhiValues, rFTmp);

      // Free temporary data
      this->shManipulator().freeTmp(rFTmp);
   }

   template <typename TSHTraits> void TorPolSHTransform<TSHTraits>::transformST2Phi(RTPScalar &rPhiValues, const typename TorPolSHTransform<TSHTraits>::SHScalarType &sComp, const typename TorPolSHTransform<TSHTraits>::SHScalarType &tComp)
   {
      if(this->shManipulator().isFirstEntry())
      {
         // add Packet size of transform
         this->shManipulator().addPacks(1);
      }

      if(this->shManipulator().atSendStage())
      {
         // add Packet size of transform
         FFTBackwardType &rMTmp = this->shManipulator().provideBTmp();

         // Perform summation over Legendre polynomials: compute S contribution
         this->legT().template multM<SetProduct>(rMTmp, sComp, &AssociatedLegendreOperator::projS2Ph);

         // Perform summation over Legendre polynomials: comput T contribution
         this->legT().template multM<AddProduct>(rMTmp, tComp, &AssociatedLegendreOperator::projT2Ph);

         // Send data and release temporary storage
         this->shManipulator().send(rMTmp);
         this->shManipulator().releaseTmp(rMTmp);
      }

      if(this->shManipulator().atRecvStage())
      {
         this->recv_ST2Phi(rPhiValues);
      }
   }

   template <typename TSHTraits> void TorPolSHTransform<TSHTraits>::transformTP2T(const RTPScalar &thetaComp, const RTPScalar &phiComp)
   {
      // add Packet size of transform
      this->shManipulator().addPacks(2);

      this->send_TP2T(thetaComp, phiComp);
   }

   template <typename TSHTraits> void TorPolSHTransform<TSHTraits>::send_TP2T(const RTPScalar &thetaComp, const RTPScalar &phiComp)
   {
      // Get temporary storage
      FFTForwardType  &rFTmpTh = this->shManipulator().provideFTmp();

      // Do the FFT transform of the real space values
      this->fft().forward(rFTmpTh, thetaComp);

      // Send data and release temporary storage
      this->shManipulator().send(rFTmpTh);
      this->shManipulator().releaseTmp(rFTmpTh);

      // Get temporary storage
      FFTForwardType  &rFTmpPh = this->shManipulator().provideFTmp();

      // Do the FFT transform of the real space values
      this->fft().forward(rFTmpPh, phiComp);

      // Send data and release temporary storage
      this->shManipulator().send(rFTmpPh);
      this->shManipulator().releaseTmp(rFTmpPh);
   }

   template <typename TSHTraits> void TorPolSHTransform<TSHTraits>::transformTP2T(typename TorPolSHTransform<TSHTraits>::SHScalarType &rTComp, const RTPScalar &thetaComp, const RTPScalar &phiComp)
   {
      if(this->shManipulator().isFirstEntry())
      {
         // add Packet size of transform
         this->shManipulator().addPacks(2);
      }

      if(this->shManipulator().atSendStage())
      {
         this->send_TP2T(thetaComp, phiComp);
      }

      if(this->shManipulator().atRecvStage())
      {
         // Get temporary storage and receive data
         FFTBackwardType &rMTmpTh = this->shManipulator().provideBTmp();
         this->shManipulator().receive(rMTmpTh);

         // Perform integration over Legendre polynomials: compute T component
         this->legT().template multM<SetProduct>(rTComp, rMTmpTh, &AssociatedLegendreOperator::intgTh2T);

         // Free temporary data
         this->shManipulator().freeTmp(rMTmpTh);

         // Get temporary storage and receive data
         FFTBackwardType &rMTmpPh = this->shManipulator().provideBTmp();
         this->shManipulator().receive(rMTmpPh);

         // Perform integration over Legendre polynomials: compute T component
         this->legT().template multM<AddProduct>(rTComp,rMTmpPh, &AssociatedLegendreOperator::intgPh2T);

         // Free temporary data
         this->shManipulator().freeTmp(rMTmpPh);

         // Set imaginary part of m=0 mode to zero
         //    Need to check if m=0 ist stored on this cpu
         if(this->hasM0())
         {
            int nR = rTComp.nR(0);
            int nL = rTComp.nL(0);
            for(int j=0; j < nR; ++j)
            {
               for(int i=0; i < nL; ++i)
               {
                  rTComp.rMShell(0)(i,j).imag() = 0.0;
               }
            }
         }
      }
   }

   template <typename TSHTraits> void TorPolSHTransform<TSHTraits>::transformSpec2GradTP(RTPScalar &rThComp, RTPScalar &rPhComp)
   {
      this->recv_Spec2GradTP(rThComp, rPhComp);
   }

   template <typename TSHTraits> void TorPolSHTransform<TSHTraits>::recv_Spec2GradTP(RTPScalar &rThComp, RTPScalar &rPhComp)
   {
      // Get temporary storage and receive data
      FFTForwardType &rFTmpTh = this->shManipulator().provideFTmp();
      rFTmpTh.doZeroPadding();
      this->shManipulator().receive(rFTmpTh);

      // Do the FFT of the complex space values
      this->fft().backward(rThComp, rFTmpTh);

      // Free temporary data
      this->shManipulator().freeTmp(rFTmpTh);

      // Get temporary storage and receive data
      FFTForwardType &rFTmpPh = this->shManipulator().provideFTmp();
      rFTmpPh.doZeroPadding();
      this->shManipulator().receive(rFTmpPh);

      // Do the FFT of the complex space values
      this->fft().backward(rPhComp, rFTmpPh);

      // Free temporary data
      this->shManipulator().freeTmp(rFTmpPh);
   }

   template <typename TSHTraits> void TorPolSHTransform<TSHTraits>::transformSpec2GradTP(RTPScalar &rThComp, RTPScalar &rPhComp, const typename TorPolSHTransform<TSHTraits>::SHScalarType &scalar)
   {
      if(this->shManipulator().isFirstEntry())
      {
         // add Packet size of transform
         this->shManipulator().addPacks(2);
      }

      if(this->shManipulator().atSendStage())
      {
         // add Packet size of transform
         FFTBackwardType &rMTmpTh = this->shManipulator().provideBTmp();

         // Perform summation over Legendre polynomials: compute S contribution
         this->legT().template multM<SetProduct>(rMTmpTh, scalar, &AssociatedLegendreOperator::proj2GradTh);

         // Send data and release temporary storage
         this->shManipulator().send(rMTmpTh);
         this->shManipulator().releaseTmp(rMTmpTh);

         // add Packet size of transform
         FFTBackwardType &rMTmpPh = this->shManipulator().provideBTmp();

         // Perform summation over Legendre polynomials: compute S contribution
         this->legT().template multM<SetProduct>(rMTmpPh, scalar, &AssociatedLegendreOperator::proj2GradPh);

         // Send data and release temporary storage
         this->shManipulator().send(rMTmpPh);
         this->shManipulator().releaseTmp(rMTmpPh);
      }

      if(this->shManipulator().atRecvStage())
      {
         this->recv_Spec2GradTP(rThComp, rPhComp);
      }
   }

}

#endif // TORPOLSHTRANSFORME_HPP
