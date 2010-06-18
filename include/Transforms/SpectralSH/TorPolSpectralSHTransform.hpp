/** \file TorPolSpectralSHTransform.hpp
 *  \brief Implementation of the Toroidal/Poloidal spectral transforms
 */

#ifndef TORPOLSPECTRALSHTRANSFORM_HPP
#define TORPOLSPECTRALSHTRANSFORM_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"
#include "GeneralScalars/RTPScalar.hpp"
#include "Transforms/SphericalHarmonics/TorPolSHTransform.hpp"
#include "Polynomials/PolyOps/SetProduct.hpp"
#include "Polynomials/PolyOps/AddProduct.hpp"

namespace EPMDynamo {

   /**
    * @brief Implements the real space to SpectralSH representation transforms
    *
    * This transforms extend the FDSH transforms to add in a radial spectral expansion.
    * The spherical harmonics expansion is dealt with by the SHTransform class. This
    * extension is responsible for doing the radial expansion on the worland polynomials.
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSpectralSHTraits Traits for radial expansion transforms
    * \tparam TSHTraits Traits for the spherical harmonics transforms
    *
    * \bug Review Implementation and traits, try to simplify flow control
    */
   template <typename TSimType, template <typename> class TSpectralSHTraits, typename TSHTraits> class TorPolSpectralSHTransform : public TorPolSHTransform<TSHTraits>
   {
      public:
         /// Typedef for the Spectral SH scalar data type
         typedef typename TSimType::ScalarType  SpectralSHScalarType;

         /// Typedef for the Spectral radial Transform data type
         typedef typename TSimType::RadialTransformType RadialTransformType;

         /// Typedef for the FDSH forward data type scalar data type
         typedef typename TSpectralSHTraits<TSimType>::FDSHForwardType  FDSHForwardType;

         /// Typedef for the FDSH backward data type
         typedef typename TSpectralSHTraits<TSimType>::FDSHBackwardType FDSHBackwardType;

         /// Typedef for the DataManipulator data type
         typedef typename TSpectralSHTraits<TSimType>::DataManipulatorType   DataManipulatorType;

         /// Typedef for the Spectral radial Transform data type
         typedef typename RadialTransformType::PolynomialBasis RadialBasisType;

         /// Typedef for the operator in the basis
         typedef typename RadialBasisType::PolynomialType  PolynomialType;

         /// Typedef for the Spectral Toroidal/Poloidal data type
         typedef typename TSpectralSHTraits<TSimType>::TorPolFieldType TorPolFieldType;

         /// Typedef for the SH Transform data type
         typedef TorPolSHTransform<TSHTraits> SHTransform;

         /**
         * @brief Constructor the SpectralSH transform and underlying SH transform
         *
         * This constructor will build the two transforms required for the Spherical
         * harmonics expansion and add to it the radial expansion: an FFT, a Legendre
         * transform and a Worland transform.
         *
         * @param pTrunc Full spectral truncation information
         */
         TorPolSpectralSHTransform(const SmartTruncation pTrunc);

         /**
         * @brief Simple empty destructor
         */
         virtual ~TorPolSpectralSHTransform() {};

         /**
          * @brief Transform RTP into SpectralSH decomposition
          *
          * @param rSSHValues Ouput full spectral decomposition coefficients
          * @param rtpValues  Input real space values
          */
         void transformRTP2Spec(SpectralSHScalarType   &rSSHValues, const RTPScalar   &rtpValues);

         /**
          * @brief Transform SpectralSH decomposition into RTP
          *
          * @param rRTPValues  Output real space values
          * @param sshValues Input full spectral decomposition coefficients
          */
         void transformSpec2RTP(RTPScalar   &rRTPValues, const SpectralSHScalarType   &sshValues);

         /**
          * @brief Transform SpectralSH decomposition into RTP gradient
          *
          * @param rGradient  Output real space gradient values
          * @param sshValues Input full spectral decomposition coefficients
          */
         void transformSpec2Grad(RTPField   &rGradient, const SpectralSHScalarType   &sshValues);

         /**
          * @brief Transform RTP into Toroidal/Poloidal decomposition
          *
          * @param rToPoField  Output Toroidal Poloida field
          * @param rtpField  Input real space values
          */
         void transformRTP2TorPol(TorPolFieldType &rToPoField, const RTPField &rtpField);

         /**
          * @brief Transform Toroidal/Poloidal decomposition into RTP
          *
          * @param rRTPField Output real space values
          * @param topoField  Input Toroidal Poloidal field
          */
         void transformTorPol2RTP(RTPField   &rRTPField, const TorPolFieldType &topoField);

         /**
          * @brief Transform Toroidal/Poloidal decomposition into RTP curl values
          *
          * @param rCurl Output real space values
          * @param topoField  Input Toroidal poloidal field
          */
         void transformTorPol2Curl(RTPField   &rCurl, const TorPolFieldType &topoField);

         /**
          * @brief Transform RTP into Non linear terms for Toroidal/Poloidal equation
          *
          * This transforms computes the following two values
          *    \f[
          *       \frac{r}{l(l+1)}\hat{r}\cdot\nabla\times\vec{N}
          *    \f]
          *    and 
          *    \f[
          *       -\frac{r}{l(l+1)}\hat{r}\cdot\nabla\times\nabla\times\vec{N}
          *    \f]
          *
          * @param rCurl  Output curl of rtp field
          * @param rCurlCurl  Output double curl of rtp field
          * @param rtpField  Input real space values
          */
         void transformRTP2TorPolNTerms(SpectralSHScalarType &rCurl, SpectralSHScalarType &rCurlCurl, const RTPField &rtpField);

         /**
          * @brief Get the radial basis
          */
         RadialBasisType&  radBasis();

         /**
          * @brief Get the radial basis
          */
         const RadialBasisType&  radBasis() const;

         /**
          * @brief Finishes initialisation of RTP Domains by filling grid values
          *
          * @param rTrunc Truncation information
          */
         void initRTPDomains(SmartTruncation &rTrunc) const;

         /**
          * @brief Get the data manipulator
          */
         DataManipulatorType&  spectralSHManipulator();
         
      protected:

         /**
          * @brief Transform Poloidal Component into R component
          *
          * @param rRComp  R Component from RTP decomposition
          */
         void transformPol2R(RTPScalar   &rRComp);

         /**
          * @brief Transform Poloidal Component into R component
          *
          * @param rRComp  R Component from RTP decomposition
          * @param polComp Poloidal Component
          */
         void transformPol2R(RTPScalar   &rRComp, const SpectralSHScalarType   &polComp);

      private:
         /**
          * @brief Radial transform object for the radial part
          */
         RadialTransformType   mSpectralTransform;

         /**
          * @brief Data manipulator
          */
         DataManipulatorType   mSpectralSHManipulator;

         /**
          * @brief Get the spectral transform
          */
         RadialTransformType&  spectralT();
   };

   template <typename TSimType, template <typename> class TSpectralSHTraits, typename TSHTraits> inline typename TorPolSpectralSHTransform<TSimType, TSpectralSHTraits, TSHTraits>::DataManipulatorType& TorPolSpectralSHTransform<TSimType, TSpectralSHTraits, TSHTraits>::spectralSHManipulator()
   {
      return this->mSpectralSHManipulator;
   }

   template <typename TSimType, template <typename> class TSpectralSHTraits, typename TSHTraits> inline typename TorPolSpectralSHTransform<TSimType, TSpectralSHTraits, TSHTraits>::RadialBasisType& TorPolSpectralSHTransform<TSimType, TSpectralSHTraits, TSHTraits>::radBasis()
   {
      return this->spectralT().basis();
   }

   template <typename TSimType, template <typename> class TSpectralSHTraits, typename TSHTraits> inline const typename TorPolSpectralSHTransform<TSimType, TSpectralSHTraits, TSHTraits>::RadialBasisType& TorPolSpectralSHTransform<TSimType, TSpectralSHTraits, TSHTraits>::radBasis() const
   {
      return this->mSpectralTransform.basis();
   }

   template <typename TSimType, template <typename> class TSpectralSHTraits, typename TSHTraits> inline typename TorPolSpectralSHTransform<TSimType, TSpectralSHTraits, TSHTraits>::RadialTransformType& TorPolSpectralSHTransform<TSimType, TSpectralSHTraits, TSHTraits>::spectralT()
   {
      return this->mSpectralTransform;
   }

   template <typename TSimType, template <typename> class TSpectralSHTraits, typename TSHTraits> TorPolSpectralSHTransform<TSimType, TSpectralSHTraits, TSHTraits>::TorPolSpectralSHTransform(const SmartTruncation pTrunc)
      : TorPolSHTransform<TSHTraits>(pTrunc), mSpectralTransform(pTrunc), mSpectralSHManipulator(pTrunc)
   {
   }

   template <typename TSimType, template <typename> class TSpectralSHTraits, typename TSHTraits> void TorPolSpectralSHTransform<TSimType, TSpectralSHTraits, TSHTraits>::transformRTP2Spec(typename TorPolSpectralSHTransform<TSimType, TSpectralSHTraits, TSHTraits>::SpectralSHScalarType &rSSHValues, const RTPScalar &rtpValues)
   {
      if(this->spectralSHManipulator().isFirstEntry())
      {
         // add Packet size of transform
         this->spectralSHManipulator().addPacks(1);
      }

      if(this->spectralSHManipulator().atSpecialStage())
      {
         // Perform the SH expansion and store in temporary FDSHForwardType
         SHTransform::transformRTP2Spec(rtpValues);
      }

      if(this->spectralSHManipulator().atSendStage())
      {
         // Get temporary storage
         FDSHForwardType  &rMTmp = this->spectralSHManipulator().provideFTmp();

         // Perform the SH expansion and store in temporary FDSHForwardType
         SHTransform::transformRTP2Spec(rMTmp, rtpValues);

         // Send data and release temporary storage
         this->spectralSHManipulator().send(rMTmp);
         this->spectralSHManipulator().releaseTmp(rMTmp);
      }

      if(this->spectralSHManipulator().atRecvStage())
      {
         // Get temporary storage and receive data
         FDSHBackwardType &rLTmp = this->spectralSHManipulator().provideBTmp();
         this->spectralSHManipulator().receive(rLTmp);

         // Perform the forward Worland transform 
         this->spectralT().template multL<SetProduct>(rSSHValues, rLTmp, &PolynomialType::intg);

         // Free temporary data
         this->spectralSHManipulator().freeTmp(rLTmp);
      }
   }

   template <typename TSimType, template <typename> class TSpectralSHTraits, typename TSHTraits> void TorPolSpectralSHTransform<TSimType, TSpectralSHTraits, TSHTraits>::transformSpec2RTP(RTPScalar &rRTPValues, const typename TorPolSpectralSHTransform<TSimType, TSpectralSHTraits, TSHTraits>::SpectralSHScalarType &sshValues)
   {
      if(this->spectralSHManipulator().isFirstEntry())
      {
         // add Packet size of transform
         this->spectralSHManipulator().addPacks(1);
      }

      if(this->spectralSHManipulator().atSendStage())
      {
         // add Packet size of transform
         FDSHBackwardType &rLTmp = this->spectralSHManipulator().provideBTmp();

         // Perform the backward Worland transform
         this->spectralT().template multL<SetProduct>(rLTmp, sshValues, &PolynomialType::proj);

         // Send data and release temporary storage
         this->spectralSHManipulator().send(rLTmp);
         this->spectralSHManipulator().releaseTmp(rLTmp);
      }

      if(this->spectralSHManipulator().atRecvStage())
      {
         // Get temporary storage and receive data
         FDSHForwardType &rMTmp = this->spectralSHManipulator().provideFTmp();
         this->spectralSHManipulator().receive(rMTmp);

         // Extract real space values from the SH expansion
         SHTransform::transformSpec2RTP(rRTPValues, rMTmp);

         // Free temporary data
         this->spectralSHManipulator().freeTmp(rMTmp);
      }

      if(this->spectralSHManipulator().atSpecialStage())
      {
         // Extract real space values from the SH expansion
         SHTransform::transformSpec2RTP(rRTPValues);
      }
   }

   template <typename TSimType, template <typename> class TSpectralSHTraits, typename TSHTraits> void TorPolSpectralSHTransform<TSimType, TSpectralSHTraits, TSHTraits>::transformSpec2Grad(RTPField &rGradient, const typename TorPolSpectralSHTransform<TSimType, TSpectralSHTraits, TSHTraits>::SpectralSHScalarType &sshValues)
   {
      if(this->spectralSHManipulator().isFirstEntry())
      {
         // add Packet size of transform
         this->spectralSHManipulator().addPacks(2);
      }

      if(this->spectralSHManipulator().atSendStage())
      {
         // Get temporary storage
         FDSHBackwardType &rLTmpR = this->spectralSHManipulator().provideBTmp();

         // Perform the backward Worland transform to obtain radial derivative
         this->spectralT().template multL<SetProduct>(rLTmpR, sshValues, &PolynomialType::dProj);

         // Transpose the FDSHForwardType values to a suitable order for the SHTransform
         this->spectralSHManipulator().send(rLTmpR);
         this->spectralSHManipulator().releaseTmp(rLTmpR);

         // Get temporary storage
         FDSHBackwardType &rLTmpTP = this->spectralSHManipulator().provideBTmp();

         // Perform the backward Worland transform to obtain value divided by the radial grid
         this->spectralT().template multL<SetProduct>(rLTmpTP, sshValues, &PolynomialType::proj2GradTP);

         // Transpose the FDSHForwardType values to a suitable order for the SHTransform
         this->spectralSHManipulator().send(rLTmpTP);
         this->spectralSHManipulator().releaseTmp(rLTmpTP);
      }

      if(this->spectralSHManipulator().atRecvStage())
      {
         // Get temporary storage and receive data
         FDSHForwardType &rMTmpR = this->spectralSHManipulator().provideFTmp();
         this->spectralSHManipulator().receive(rMTmpR);

         // Extract real space values from the SH expansion
         SHTransform::transformSpec2RTP(rGradient.rR(), rMTmpR);

         // Free temporary data
         this->spectralSHManipulator().freeTmp(rMTmpR);

         // Get temporary storage and receive data
         FDSHForwardType &rMTmpTP = this->spectralSHManipulator().provideFTmp();
         this->spectralSHManipulator().receive(rMTmpTP);

         // Extract real space values from the SH expansion
         SHTransform::transformSpec2GradTP(rGradient.rTheta(), rGradient.rPhi(), rMTmpTP);

         // Free temporary data
         this->spectralSHManipulator().freeTmp(rMTmpTP);
      }

      if(this->spectralSHManipulator().atSpecialStage())
      {
         SHTransform::transformSpec2RTP(rGradient.rR());

         SHTransform::transformSpec2GradTP(rGradient.rTheta(), rGradient.rPhi());
      }
   }

   template <typename TSimType, template <typename> class TSpectralSHTraits, typename TSHTraits> void TorPolSpectralSHTransform<TSimType, TSpectralSHTraits, TSHTraits>::transformRTP2TorPol(TorPolFieldType &rToPoField, const RTPField &rtpField)
   {
      if(this->spectralSHManipulator().isFirstEntry())
      {
         // add Packet size of transform
         this->spectralSHManipulator().addPacks(2);
      }

      if(this->spectralSHManipulator().atSpecialStage())
      {
         // Compute Q component out of r component
         SHTransform::transformRTP2Spec(rtpField.r());

         // Compute T component from the Theta and phi components
         this->transformTP2T(rtpField.theta(), rtpField.phi());
      }

      if(this->spectralSHManipulator().atSendStage())
      {
         // Get temporary storage
         FDSHForwardType &rMTmpQ = this->spectralSHManipulator().provideFTmp();

         // Compute Q component out of r component
         SHTransform::transformRTP2Spec(rMTmpQ, rtpField.r());

         // Transpose the FDSHForwardType values to a suitable order for radial transform
         this->spectralSHManipulator().send(rMTmpQ);
         this->spectralSHManipulator().releaseTmp(rMTmpQ);

         // Get temporary storage
         FDSHForwardType &rMTmpT = this->spectralSHManipulator().provideFTmp();

         // Compute T component from the Theta and phi components
         this->transformTP2T(rMTmpT, rtpField.theta(), rtpField.phi());

         // Transpose the FDSHForwardType values to a suitable order for radial transform
         this->spectralSHManipulator().send(rMTmpT);
         this->spectralSHManipulator().releaseTmp(rMTmpT);
      }

      if(this->spectralSHManipulator().atRecvStage())
      {
         // Get temporary storage and receive data
         FDSHBackwardType &rLTmpQ = this->spectralSHManipulator().provideBTmp();
         this->spectralSHManipulator().receive(rLTmpQ);

         // Perform the forward Worland transform for Poloidal component
         this->spectralT().template multL<SetProduct>(rToPoField.rPol(), rLTmpQ, &PolynomialType::intgQ2Pol);

         // Free temporary data
         this->spectralSHManipulator().freeTmp(rLTmpQ);

         // Get temporary storage and receive data
         FDSHBackwardType &rLTmpT = this->spectralSHManipulator().provideBTmp();
         this->spectralSHManipulator().receive(rLTmpT);

         // Perform the forward Worland transform for Toroidal component
         this->spectralT().template multL<SetProduct>(rToPoField.rTor(), rLTmpT, &PolynomialType::intgT2Tor);

         // Free temporary data
         this->spectralSHManipulator().freeTmp(rLTmpT);
      }
   }

   template <typename TSimType, template <typename> class TSpectralSHTraits, typename TSHTraits> void TorPolSpectralSHTransform<TSimType, TSpectralSHTraits, TSHTraits>::transformRTP2TorPolNTerms(typename TorPolSpectralSHTransform<TSimType, TSpectralSHTraits, TSHTraits>::SpectralSHScalarType &rCurl, typename TorPolSpectralSHTransform<TSimType, TSpectralSHTraits, TSHTraits>::SpectralSHScalarType &rCurlCurl, const RTPField &rtpField)
   {
      if(this->spectralSHManipulator().isFirstEntry())
      {
         // add Packet size of transform
         this->spectralSHManipulator().addPacks(3);
      }

      if(this->spectralSHManipulator().atSpecialStage())
      {
         // Get the Q component from radial component
         SHTransform::transformRTP2Spec(rtpField.r());

         // Compute S and T components: Theta part
         this->transformTP2STPartT(rtpField.theta());

         // Compute S and T components: Phi part
         this->transformTP2STPartP(rtpField.phi());
      }

      if(this->spectralSHManipulator().atSendStage())
      {
         // Get temporary storage
         FDSHForwardType &rMTmpQ = this->spectralSHManipulator().provideFTmp();

         // Get the Q component from radial component
         SHTransform::transformRTP2Spec(rMTmpQ, rtpField.r());

         // Transpose the FDSHForwardType values to a suitable order for radial transform
         this->spectralSHManipulator().send(rMTmpQ);
         this->spectralSHManipulator().releaseTmp(rMTmpQ);

         // Get temporary storage
         FDSHForwardType &rMTmpS = this->spectralSHManipulator().provideFTmp();
         FDSHForwardType &rMTmpT = this->spectralSHManipulator().provideFTmp();

         // Compute S and T components: Theta part
         this->transformTP2STPartT(rMTmpS, rMTmpT, rtpField.theta());

         // Compute S and T components: Phi part
         this->transformTP2STPartP(rMTmpS, rMTmpT, rtpField.phi());

         // Transpose the FDSHForwardType values to a suitable order for radial transform
         this->spectralSHManipulator().send(rMTmpS);
         this->spectralSHManipulator().releaseTmp(rMTmpS);

         // Transpose the FDSHForwardType values to a suitable order for radial transform
         this->spectralSHManipulator().send(rMTmpT);
         this->spectralSHManipulator().releaseTmp(rMTmpT);
      }

      if(this->spectralSHManipulator().atRecvStage())
      {
         // Get temporary storage and receive data
         FDSHBackwardType &rLTmpQ = this->spectralSHManipulator().provideBTmp();
         this->spectralSHManipulator().receive(rLTmpQ);

         // Perform the forward Worland transform  from Q component
         this->spectralT().template multL<SetProduct>(rCurlCurl, rLTmpQ, &PolynomialType::intgQ2CurlCurlProj);
         
         // Free temporary storage
         this->spectralSHManipulator().freeTmp(rLTmpQ);

         // Get temporary storage and receive data
         FDSHBackwardType &rLTmpS = this->spectralSHManipulator().provideBTmp();
         this->spectralSHManipulator().receive(rLTmpS);

         // Perform the forward Worland transform  from S component
         this->spectralT().template multL<AddProduct>(rCurlCurl, rLTmpS, &PolynomialType::intgS2CurlCurlProj);

         // Free temporary storage
         this->spectralSHManipulator().freeTmp(rLTmpS);

         // Get temporary storage and receive data
         FDSHBackwardType &rLTmpT = this->spectralSHManipulator().provideBTmp();
         this->spectralSHManipulator().receive(rLTmpT);

         // Perform the forward Worland transform for Curl value from T component
         this->spectralT().template multL<SetProduct>(rCurl, rLTmpT, &PolynomialType::intgT2CurlProj);

         // Free temporary storage
         this->spectralSHManipulator().freeTmp(rLTmpT);
      }
   }

   template <typename TSimType, template <typename> class TSpectralSHTraits, typename TSHTraits> void TorPolSpectralSHTransform<TSimType, TSpectralSHTraits, TSHTraits>::transformTorPol2RTP(RTPField &rRTPField, const TorPolFieldType &topoField)
   {
      if(this->spectralSHManipulator().isFirstEntry())
      {
         // add Packet size of transform
         this->spectralSHManipulator().addPacks(3);
      }

      if(this->spectralSHManipulator().atSendStage())
      {
         // Compute the R component from the poloidal component
         this->transformPol2R(rRTPField.rR(), topoField.pol());

         // Get temporary storage
         FDSHBackwardType &rLTmpS = this->spectralSHManipulator().provideBTmp();

         // Perform the backward Worland transform for Poloidal component: S 
         this->spectralT().template multL<SetProduct>(rLTmpS, topoField.pol(), &PolynomialType::projPol2S);

         // Transpose the FDSHForwardType values to a suitable order for the SHTransform
         this->spectralSHManipulator().send(rLTmpS);
         this->spectralSHManipulator().releaseTmp(rLTmpS);

         // Get temporary storage
         FDSHBackwardType &rLTmpT = this->spectralSHManipulator().provideBTmp();

         // Perform the backward Worland transform for Toroidal component
         this->spectralT().template multL<SetProduct>(rLTmpT, topoField.tor(), &PolynomialType::projTor2T);

         // Transpose the FDSHForwardType values to a suitable order for the SHTransform
         this->spectralSHManipulator().send(rLTmpT);
         this->spectralSHManipulator().releaseTmp(rLTmpT);
      }

      if(this->spectralSHManipulator().atRecvStage())
      {
         // Compute the R component from the poloidal component
         this->transformPol2R(rRTPField.rR(), topoField.pol());

         // Get temporary storage and receive data
         FDSHForwardType &rMTmpS = this->spectralSHManipulator().provideFTmp();
         this->spectralSHManipulator().receive(rMTmpS);
         FDSHForwardType &rMTmpT = this->spectralSHManipulator().provideFTmp();
         this->spectralSHManipulator().receive(rMTmpT);

         // Compute Theta component
         this->transformST2Theta(rRTPField.rTheta(), rMTmpS, rMTmpT);

         // Compute Phi component
         this->transformST2Phi(rRTPField.rPhi(), rMTmpS, rMTmpT);

         // Free temporary storage 
         this->spectralSHManipulator().freeTmp(rMTmpS);
         this->spectralSHManipulator().freeTmp(rMTmpT);
      }

      if(this->spectralSHManipulator().atSpecialStage())
      {
         // Compute the R component from the poloidal component
         this->transformPol2R(rRTPField.rR());

         // Compute Theta component
         this->transformST2Theta(rRTPField.rTheta());

         // Compute Phi component
         this->transformST2Phi(rRTPField.rPhi());
      }
   }

   template <typename TSimType, template <typename> class TSpectralSHTraits, typename TSHTraits> void TorPolSpectralSHTransform<TSimType, TSpectralSHTraits, TSHTraits>::transformTorPol2Curl(RTPField &rCurl, const TorPolFieldType &topoField)
   {
      if(this->spectralSHManipulator().isFirstEntry())
      {
         // add Packet size of transform
         this->spectralSHManipulator().addPacks(3);
      }

      if(this->spectralSHManipulator().atSendStage())
      {
         // Get temporary storage
         FDSHBackwardType &rLTmpCurlQ = this->spectralSHManipulator().provideBTmp();

         // Perform the backward Worland transform for S component
         this->spectralT().template multL<SetProduct>(rLTmpCurlQ, topoField.tor(), &PolynomialType::projTor2CurlQ);

         // Transpose the FDSHForwardType values to a suitable order for the SHTransform
         this->spectralSHManipulator().send(rLTmpCurlQ);
         this->spectralSHManipulator().releaseTmp(rLTmpCurlQ);

         // Get temporary storage
         FDSHBackwardType &rLTmpCurlS = this->spectralSHManipulator().provideBTmp();

         // Perform the backward Worland transform for S component
         this->spectralT().template multL<SetProduct>(rLTmpCurlS, topoField.tor(), &PolynomialType::projTor2CurlS);

         // Transpose the FDSHForwardType values to a suitable order for the SHTransform
         this->spectralSHManipulator().send(rLTmpCurlS);
         this->spectralSHManipulator().releaseTmp(rLTmpCurlS);

         // Get temporary storage
         FDSHBackwardType &rLTmpCurlT = this->spectralSHManipulator().provideBTmp();

         // Perform the backward Worland transform for T component
         this->spectralT().template multL<SetProduct>(rLTmpCurlT, topoField.pol(), &PolynomialType::projPol2CurlT);

         // Transpose the FDSHForwardType values to a suitable order for the SHTransform
         this->spectralSHManipulator().send(rLTmpCurlT);
         this->spectralSHManipulator().releaseTmp(rLTmpCurlT);
      }

      if(this->spectralSHManipulator().atRecvStage())
      {
         // Get temporary storage and receive data
         FDSHForwardType &rMTmpCurlQ = this->spectralSHManipulator().provideFTmp();
         this->spectralSHManipulator().receive(rMTmpCurlQ);

         // Extract real space radial component from the SH expansion
         SHTransform::transformSpec2RTP(rCurl.rR(), rMTmpCurlQ);

         // Free temporary data
         this->spectralSHManipulator().freeTmp(rMTmpCurlQ);

         // Get temporary storage and receive data
         FDSHForwardType &rMTmpCurlS = this->spectralSHManipulator().provideFTmp();
         this->spectralSHManipulator().receive(rMTmpCurlS);
         FDSHForwardType &rMTmpCurlT = this->spectralSHManipulator().provideFTmp();
         this->spectralSHManipulator().receive(rMTmpCurlT);

         // Compute Theta component
         this->transformST2Theta(rCurl.rTheta(), rMTmpCurlS, rMTmpCurlT);

         // Compute Phi component
         this->transformST2Phi(rCurl.rPhi(), rMTmpCurlS, rMTmpCurlT);

         // Free temporary data
         this->spectralSHManipulator().freeTmp(rMTmpCurlS);
         this->spectralSHManipulator().freeTmp(rMTmpCurlT);
      }

      if(this->spectralSHManipulator().atSpecialStage())
      {
         // Extract real space radial component from the SH expansion
         SHTransform::transformSpec2RTP(rCurl.rR());

         // Compute Theta component
         this->transformST2Theta(rCurl.rTheta());

         // Compute Phi component
         this->transformST2Phi(rCurl.rPhi());
      }
   }

   template <typename TSimType, template <typename> class TSpectralSHTraits, typename TSHTraits> void TorPolSpectralSHTransform<TSimType, TSpectralSHTraits, TSHTraits>::initRTPDomains(SmartTruncation &rTrunc) const
   {
      rTrunc->setRadialGrid(this->radBasis().at(0).pGrid());
      rTrunc->setThetaGrid(this->hozBasis().at(0).pTheta(), this->hozBasis().at(0).pCosTheta(), this->hozBasis().at(0).pSinTheta());
   }

   template <typename TSimType, template <typename> class TSpectralSHTraits, typename TSHTraits> void TorPolSpectralSHTransform<TSimType, TSpectralSHTraits, TSHTraits>::transformPol2R(RTPScalar &rRComp)
   {
      // Extract real space values from the SH expansion
      SHTransform::transformSpec2RTP(rRComp);
   }

   template <typename TSimType, template <typename> class TSpectralSHTraits, typename TSHTraits> void TorPolSpectralSHTransform<TSimType, TSpectralSHTraits, TSHTraits>::transformPol2R(RTPScalar &rRComp, const typename TorPolSpectralSHTransform<TSimType, TSpectralSHTraits, TSHTraits>::SpectralSHScalarType &polComp)
   {
      if(this->spectralSHManipulator().atSendStage())
      {
         // Get temporary storage
         FDSHBackwardType &rLTmp = this->spectralSHManipulator().provideBTmp();

         // Perform the backward Worland transform for Poloidal component: Q part
         this->spectralT().template multL<SetProduct>(rLTmp, polComp, &PolynomialType::projPol2Q);

         // Transpose the FDSHForwardType values to a suitable order for the SHTransform
         this->spectralSHManipulator().send(rLTmp);
         this->spectralSHManipulator().releaseTmp(rLTmp);
      } 

      if(this->spectralSHManipulator().atRecvStage())
      {
         // Get temporary storage and receive data
         FDSHForwardType &rMTmp = this->spectralSHManipulator().provideFTmp();
         this->spectralSHManipulator().receive(rMTmp);

         // Extract real space values from the SH expansion
         SHTransform::transformSpec2RTP(rRComp, rMTmp);

         // free temporary storage
         this->spectralSHManipulator().freeTmp(rMTmp);
      }
   }

}

#endif // TORPOLSPECTRALSHTRANSFORM_HPP
