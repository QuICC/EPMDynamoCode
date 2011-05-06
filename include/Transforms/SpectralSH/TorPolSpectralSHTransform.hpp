/** \file TorPolSpectralSHTransform.hpp
 *  \brief Implementation of the Toroidal/Poloidal spectral transforms
 */

#ifndef TORPOLSPECTRALSHTRANSFORM_HPP
#define TORPOLSPECTRALSHTRANSFORM_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

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
    * \tparam TSpectralSHTraits Traits for radial expansion transforms
    * \tparam TSHTraits Traits for the spherical harmonics transforms
    */
   template <typename TSpectralSHTraits, typename TSHTraits> class TorPolSpectralSHTransform : public TorPolSHTransform<TSHTraits>
   {
      public:
         /// Typedef for the Spectral SH scalar data type
         typedef SimulationConfig::NumericalScheme::ScalarType  SpectralSHScalarType;

         /// Typedef for the Spectral radial Transform data type
         typedef SimulationConfig::NumericalScheme::RadialTransformType RadialTransformType;

         /// Typedef for the FDSH forward data type scalar data type
         typedef typename TSpectralSHTraits::FDSHForwardType  FDSHForwardType;

         /// Typedef for the FDSH backward data type
         typedef typename TSpectralSHTraits::FDSHBackwardType FDSHBackwardType;

         /// Typedef for the DataManipulator data type
         typedef typename TSpectralSHTraits::DataManipulatorType   DataManipulatorType;

         /// Typedef for the Spectral radial Transform data type
         typedef typename RadialTransformType::PolynomialBasis RadialBasisType;

         /// Typedef for the operator in the basis
         typedef typename RadialBasisType::PolynomialType  PolynomialType;

         /// Typedef for the Spectral Toroidal/Poloidal data type
         typedef typename TSpectralSHTraits::TorPolFieldType TorPolFieldType;

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
         DataManipulatorType&  sshManipulator();
         
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

   template <typename TSpectralSHTraits, typename TSHTraits> inline typename TorPolSpectralSHTransform<TSpectralSHTraits, TSHTraits>::DataManipulatorType& TorPolSpectralSHTransform<TSpectralSHTraits, TSHTraits>::sshManipulator()
   {
      return this->mSpectralSHManipulator;
   }

   template <typename TSpectralSHTraits, typename TSHTraits> inline typename TorPolSpectralSHTransform<TSpectralSHTraits, TSHTraits>::RadialBasisType& TorPolSpectralSHTransform<TSpectralSHTraits, TSHTraits>::radBasis()
   {
      return this->spectralT().basis();
   }

   template <typename TSpectralSHTraits, typename TSHTraits> inline const typename TorPolSpectralSHTransform<TSpectralSHTraits, TSHTraits>::RadialBasisType& TorPolSpectralSHTransform<TSpectralSHTraits, TSHTraits>::radBasis() const
   {
      return this->mSpectralTransform.basis();
   }

   template <typename TSpectralSHTraits, typename TSHTraits> inline typename TorPolSpectralSHTransform<TSpectralSHTraits, TSHTraits>::RadialTransformType& TorPolSpectralSHTransform<TSpectralSHTraits, TSHTraits>::spectralT()
   {
      return this->mSpectralTransform;
   }

   template <typename TSpectralSHTraits, typename TSHTraits> TorPolSpectralSHTransform<TSpectralSHTraits, TSHTraits>::TorPolSpectralSHTransform(const SmartTruncation pTrunc)
      : TorPolSHTransform<TSHTraits>(pTrunc), mSpectralTransform(pTrunc), mSpectralSHManipulator(pTrunc)
   {
   }

   template <typename TSpectralSHTraits, typename TSHTraits> void TorPolSpectralSHTransform<TSpectralSHTraits, TSHTraits>::transformRTP2Spec(typename TorPolSpectralSHTransform<TSpectralSHTraits, TSHTraits>::SpectralSHScalarType &rSSHValues, const RTPScalar &rtpValues)
   {
      if(this->sshManipulator().isFirstEntry())
      {
         // add Packet size of transform
         this->sshManipulator().addPacks(1);
      }

      if(this->sshManipulator().atInterStage())
      {
         // Perform the SH expansion and store in temporary FDSHForwardType
         SHTransform::transformRTP2Spec(rtpValues);
      }

      if(this->sshManipulator().atSendStage())
      {
         // Get temporary storage
         FDSHForwardType  &rMTmp = this->sshManipulator().provideFTmp();

         // Perform the SH expansion and store in temporary FDSHForwardType
         SHTransform::transformRTP2Spec(rMTmp, rtpValues);

         // Send data and release temporary storage
         this->sshManipulator().send(rMTmp);
         this->sshManipulator().releaseTmp(rMTmp);
      }

      if(this->sshManipulator().atRecvStage())
      {
         // Get temporary storage and receive data
         FDSHBackwardType &rLTmp = this->sshManipulator().provideBTmp();
         this->sshManipulator().receive(rLTmp);

         // Perform the forward Worland transform 
         this->spectralT().template multL<SetProduct>(rSSHValues, rLTmp, &PolynomialType::intg);

         // Free temporary data
         this->sshManipulator().freeTmp(rLTmp);
      }
   }

   template <typename TSpectralSHTraits, typename TSHTraits> void TorPolSpectralSHTransform<TSpectralSHTraits, TSHTraits>::transformSpec2RTP(RTPScalar &rRTPValues, const typename TorPolSpectralSHTransform<TSpectralSHTraits, TSHTraits>::SpectralSHScalarType &sshValues)
   {
      if(this->sshManipulator().isFirstEntry())
      {
         // add Packet size of transform
         this->sshManipulator().addPacks(1);
      }

      if(this->sshManipulator().atSendStage())
      {
         // add Packet size of transform
         FDSHBackwardType &rLTmp = this->sshManipulator().provideBTmp();

         // Perform the backward Worland transform
         this->spectralT().template multL<SetProduct>(rLTmp, sshValues, &PolynomialType::proj);

         // Send data and release temporary storage
         this->sshManipulator().send(rLTmp);
         this->sshManipulator().releaseTmp(rLTmp);
      }

      if(this->sshManipulator().atRecvStage())
      {
         // Get temporary storage and receive data
         FDSHForwardType &rMTmp = this->sshManipulator().provideFTmp();
         this->sshManipulator().receive(rMTmp);

         // Extract real space values from the SH expansion
         SHTransform::transformSpec2RTP(rRTPValues, rMTmp);

         // Free temporary data
         this->sshManipulator().freeTmp(rMTmp);
      }

      if(this->sshManipulator().atInterStage())
      {
         // Extract real space values from the SH expansion
         SHTransform::transformSpec2RTP(rRTPValues);
      }
   }

   template <typename TSpectralSHTraits, typename TSHTraits> void TorPolSpectralSHTransform<TSpectralSHTraits, TSHTraits>::transformSpec2Grad(RTPField &rGradient, const typename TorPolSpectralSHTransform<TSpectralSHTraits, TSHTraits>::SpectralSHScalarType &sshValues)
   {
      if(this->sshManipulator().isFirstEntry())
      {
         // add Packet size of transform
         this->sshManipulator().addPacks(2);
      }

      if(this->sshManipulator().atSendStage())
      {
         // Get temporary storage
         FDSHBackwardType &rLTmpR = this->sshManipulator().provideBTmp();

         // Perform the backward Worland transform to obtain radial derivative
         this->spectralT().template multL<SetProduct>(rLTmpR, sshValues, &PolynomialType::dProj);

         // Transpose the FDSHForwardType values to a suitable order for the SHTransform
         this->sshManipulator().send(rLTmpR);
         this->sshManipulator().releaseTmp(rLTmpR);

         // Get temporary storage
         FDSHBackwardType &rLTmpTP = this->sshManipulator().provideBTmp();

         // Perform the backward Worland transform to obtain value divided by the radial grid
         this->spectralT().template multL<SetProduct>(rLTmpTP, sshValues, &PolynomialType::proj2GradTP);

         // Transpose the FDSHForwardType values to a suitable order for the SHTransform
         this->sshManipulator().send(rLTmpTP);
         this->sshManipulator().releaseTmp(rLTmpTP);
      }

      if(this->sshManipulator().atRecvStage())
      {
         // Get temporary storage and receive data
         FDSHForwardType &rMTmpR = this->sshManipulator().provideFTmp();
         this->sshManipulator().receive(rMTmpR);

         // Extract real space values from the SH expansion
         SHTransform::transformSpec2RTP(rGradient.rR(), rMTmpR);

         // Free temporary data
         this->sshManipulator().freeTmp(rMTmpR);

         // Get temporary storage and receive data
         FDSHForwardType &rMTmpTP = this->sshManipulator().provideFTmp();
         this->sshManipulator().receive(rMTmpTP);

         // Extract real space values from the SH expansion
         SHTransform::transformSpec2GradTP(rGradient.rTheta(), rGradient.rPhi(), rMTmpTP);

         // Free temporary data
         this->sshManipulator().freeTmp(rMTmpTP);
      }

      if(this->sshManipulator().atInterStage())
      {
         SHTransform::transformSpec2RTP(rGradient.rR());

         SHTransform::transformSpec2GradTP(rGradient.rTheta(), rGradient.rPhi());
      }
   }

   template <typename TSpectralSHTraits, typename TSHTraits> void TorPolSpectralSHTransform<TSpectralSHTraits, TSHTraits>::transformRTP2TorPol(TorPolFieldType &rToPoField, const RTPField &rtpField)
   {
      if(this->sshManipulator().isFirstEntry())
      {
         // add Packet size of transform
         this->sshManipulator().addPacks(2);
      }

      if(this->sshManipulator().atInterStage())
      {
         // Compute Q component out of r component
         SHTransform::transformRTP2Spec(rtpField.r());

         // Compute T component from the Theta and phi components
         this->transformTP2T(rtpField.theta(), rtpField.phi());
      }

      if(this->sshManipulator().atSendStage())
      {
         // Get temporary storage
         FDSHForwardType &rMTmpQ = this->sshManipulator().provideFTmp();

         // Compute Q component out of r component
         SHTransform::transformRTP2Spec(rMTmpQ, rtpField.r());

         // Transpose the FDSHForwardType values to a suitable order for radial transform
         this->sshManipulator().send(rMTmpQ);
         this->sshManipulator().releaseTmp(rMTmpQ);

         // Get temporary storage
         FDSHForwardType &rMTmpT = this->sshManipulator().provideFTmp();

         // Compute T component from the Theta and phi components
         this->transformTP2T(rMTmpT, rtpField.theta(), rtpField.phi());

         // Transpose the FDSHForwardType values to a suitable order for radial transform
         this->sshManipulator().send(rMTmpT);
         this->sshManipulator().releaseTmp(rMTmpT);
      }

      if(this->sshManipulator().atRecvStage())
      {
         // Get temporary storage and receive data
         FDSHBackwardType &rLTmpQ = this->sshManipulator().provideBTmp();
         this->sshManipulator().receive(rLTmpQ);

         // Perform the forward Worland transform for Poloidal component
         this->spectralT().template multL<SetProduct>(rToPoField.rPol(), rLTmpQ, &PolynomialType::intgQ2Pol);

         // Free temporary data
         this->sshManipulator().freeTmp(rLTmpQ);

         // Get temporary storage and receive data
         FDSHBackwardType &rLTmpT = this->sshManipulator().provideBTmp();
         this->sshManipulator().receive(rLTmpT);

         // Perform the forward Worland transform for Toroidal component
         this->spectralT().template multL<SetProduct>(rToPoField.rTor(), rLTmpT, &PolynomialType::intgT2Tor);

         // Free temporary data
         this->sshManipulator().freeTmp(rLTmpT);
      }
   }

   template <typename TSpectralSHTraits, typename TSHTraits> void TorPolSpectralSHTransform<TSpectralSHTraits, TSHTraits>::transformRTP2TorPolNTerms(typename TorPolSpectralSHTransform<TSpectralSHTraits, TSHTraits>::SpectralSHScalarType &rCurl, typename TorPolSpectralSHTransform<TSpectralSHTraits, TSHTraits>::SpectralSHScalarType &rCurlCurl, const RTPField &rtpField)
   {
      if(this->sshManipulator().isFirstEntry())
      {
         // add Packet size of transform
         this->sshManipulator().addPacks(3);
      }

      if(this->sshManipulator().atInterStage())
      {
         // Get the Q component from radial component
         SHTransform::transformRTP2Spec(rtpField.r());

         // Compute S and T components: Theta part
         this->transformTP2STPartT(rtpField.theta());

         // Compute S and T components: Phi part
         this->transformTP2STPartP(rtpField.phi());
      }

      if(this->sshManipulator().atSendStage())
      {
         // Get temporary storage
         FDSHForwardType &rMTmpQ = this->sshManipulator().provideFTmp();

         // Get the Q component from radial component
         SHTransform::transformRTP2Spec(rMTmpQ, rtpField.r());

         // Transpose the FDSHForwardType values to a suitable order for radial transform
         this->sshManipulator().send(rMTmpQ);
         this->sshManipulator().releaseTmp(rMTmpQ);

         // Get temporary storage
         FDSHForwardType &rMTmpS = this->sshManipulator().provideFTmp();
         FDSHForwardType &rMTmpT = this->sshManipulator().provideFTmp();

         // Compute S and T components: Theta part
         this->transformTP2STPartT(rMTmpS, rMTmpT, rtpField.theta());

         // Compute S and T components: Phi part
         this->transformTP2STPartP(rMTmpS, rMTmpT, rtpField.phi());

         // Transpose the FDSHForwardType values to a suitable order for radial transform
         this->sshManipulator().send(rMTmpS);
         this->sshManipulator().releaseTmp(rMTmpS);

         // Transpose the FDSHForwardType values to a suitable order for radial transform
         this->sshManipulator().send(rMTmpT);
         this->sshManipulator().releaseTmp(rMTmpT);
      }

      if(this->sshManipulator().atRecvStage())
      {
         // Get temporary storage and receive data
         FDSHBackwardType &rLTmpQ = this->sshManipulator().provideBTmp();
         this->sshManipulator().receive(rLTmpQ);

         // Perform the forward Worland transform  from Q component
         this->spectralT().template multL<SetProduct>(rCurlCurl, rLTmpQ, &PolynomialType::intgQ2CurlCurlProj);
         
         // Free temporary storage
         this->sshManipulator().freeTmp(rLTmpQ);

         // Get temporary storage and receive data
         FDSHBackwardType &rLTmpS = this->sshManipulator().provideBTmp();
         this->sshManipulator().receive(rLTmpS);

         // Perform the forward Worland transform  from S component
         this->spectralT().template multL<AddProduct>(rCurlCurl, rLTmpS, &PolynomialType::intgS2CurlCurlProj);

         // Free temporary storage
         this->sshManipulator().freeTmp(rLTmpS);

         // Get temporary storage and receive data
         FDSHBackwardType &rLTmpT = this->sshManipulator().provideBTmp();
         this->sshManipulator().receive(rLTmpT);

         // Perform the forward Worland transform for Curl value from T component
         this->spectralT().template multL<SetProduct>(rCurl, rLTmpT, &PolynomialType::intgT2CurlProj);

         // Free temporary storage
         this->sshManipulator().freeTmp(rLTmpT);
      }
   }

   template <typename TSpectralSHTraits, typename TSHTraits> void TorPolSpectralSHTransform<TSpectralSHTraits, TSHTraits>::transformTorPol2RTP(RTPField &rRTPField, const TorPolFieldType &topoField)
   {
      if(this->sshManipulator().isFirstEntry())
      {
         // add Packet size of transform
         this->sshManipulator().addPacks(3);
      }

      if(this->sshManipulator().atSendStage())
      {
         // Compute the R component from the poloidal component
         this->transformPol2R(rRTPField.rR(), topoField.pol());

         // Get temporary storage
         FDSHBackwardType &rLTmpS = this->sshManipulator().provideBTmp();

         // Perform the backward Worland transform for Poloidal component: S 
         this->spectralT().template multL<SetProduct>(rLTmpS, topoField.pol(), &PolynomialType::projPol2S);

         // Transpose the FDSHForwardType values to a suitable order for the SHTransform
         this->sshManipulator().send(rLTmpS);
         this->sshManipulator().releaseTmp(rLTmpS);

         // Get temporary storage
         FDSHBackwardType &rLTmpT = this->sshManipulator().provideBTmp();

         // Perform the backward Worland transform for Toroidal component
         this->spectralT().template multL<SetProduct>(rLTmpT, topoField.tor(), &PolynomialType::projTor2T);

         // Transpose the FDSHForwardType values to a suitable order for the SHTransform
         this->sshManipulator().send(rLTmpT);
         this->sshManipulator().releaseTmp(rLTmpT);
      }

      if(this->sshManipulator().atRecvStage())
      {
         // Compute the R component from the poloidal component
         this->transformPol2R(rRTPField.rR(), topoField.pol());

         // Get temporary storage and receive data
         FDSHForwardType &rMTmpS = this->sshManipulator().provideFTmp();
         this->sshManipulator().receive(rMTmpS);
         FDSHForwardType &rMTmpT = this->sshManipulator().provideFTmp();
         this->sshManipulator().receive(rMTmpT);

         // Compute Theta component
         this->transformST2Theta(rRTPField.rTheta(), rMTmpS, rMTmpT);

         // Compute Phi component
         this->transformST2Phi(rRTPField.rPhi(), rMTmpS, rMTmpT);

         // Free temporary storage 
         this->sshManipulator().freeTmp(rMTmpS);
         this->sshManipulator().freeTmp(rMTmpT);
      }

      if(this->sshManipulator().atInterStage())
      {
         // Compute the R component from the poloidal component
         this->transformPol2R(rRTPField.rR());

         // Compute Theta component
         this->transformST2Theta(rRTPField.rTheta());

         // Compute Phi component
         this->transformST2Phi(rRTPField.rPhi());
      }
   }

   template <typename TSpectralSHTraits, typename TSHTraits> void TorPolSpectralSHTransform<TSpectralSHTraits, TSHTraits>::transformTorPol2Curl(RTPField &rCurl, const TorPolFieldType &topoField)
   {
      if(this->sshManipulator().isFirstEntry())
      {
         // add Packet size of transform
         this->sshManipulator().addPacks(3);
      }

      if(this->sshManipulator().atSendStage())
      {
         // Get temporary storage
         FDSHBackwardType &rLTmpCurlQ = this->sshManipulator().provideBTmp();

         // Perform the backward Worland transform for S component
         this->spectralT().template multL<SetProduct>(rLTmpCurlQ, topoField.tor(), &PolynomialType::projTor2CurlQ);

         // Transpose the FDSHForwardType values to a suitable order for the SHTransform
         this->sshManipulator().send(rLTmpCurlQ);
         this->sshManipulator().releaseTmp(rLTmpCurlQ);

         // Get temporary storage
         FDSHBackwardType &rLTmpCurlS = this->sshManipulator().provideBTmp();

         // Perform the backward Worland transform for S component
         this->spectralT().template multL<SetProduct>(rLTmpCurlS, topoField.tor(), &PolynomialType::projTor2CurlS);

         // Transpose the FDSHForwardType values to a suitable order for the SHTransform
         this->sshManipulator().send(rLTmpCurlS);
         this->sshManipulator().releaseTmp(rLTmpCurlS);

         // Get temporary storage
         FDSHBackwardType &rLTmpCurlT = this->sshManipulator().provideBTmp();

         // Perform the backward Worland transform for T component
         this->spectralT().template multL<SetProduct>(rLTmpCurlT, topoField.pol(), &PolynomialType::projPol2CurlT);

         // Transpose the FDSHForwardType values to a suitable order for the SHTransform
         this->sshManipulator().send(rLTmpCurlT);
         this->sshManipulator().releaseTmp(rLTmpCurlT);
      }

      if(this->sshManipulator().atRecvStage())
      {
         // Get temporary storage and receive data
         FDSHForwardType &rMTmpCurlQ = this->sshManipulator().provideFTmp();
         this->sshManipulator().receive(rMTmpCurlQ);

         // Extract real space radial component from the SH expansion
         SHTransform::transformSpec2RTP(rCurl.rR(), rMTmpCurlQ);

         // Free temporary data
         this->sshManipulator().freeTmp(rMTmpCurlQ);

         // Get temporary storage and receive data
         FDSHForwardType &rMTmpCurlS = this->sshManipulator().provideFTmp();
         this->sshManipulator().receive(rMTmpCurlS);
         FDSHForwardType &rMTmpCurlT = this->sshManipulator().provideFTmp();
         this->sshManipulator().receive(rMTmpCurlT);

         // Compute Theta component
         this->transformST2Theta(rCurl.rTheta(), rMTmpCurlS, rMTmpCurlT);

         // Compute Phi component
         this->transformST2Phi(rCurl.rPhi(), rMTmpCurlS, rMTmpCurlT);

         // Free temporary data
         this->sshManipulator().freeTmp(rMTmpCurlS);
         this->sshManipulator().freeTmp(rMTmpCurlT);
      }

      if(this->sshManipulator().atInterStage())
      {
         // Extract real space radial component from the SH expansion
         SHTransform::transformSpec2RTP(rCurl.rR());

         // Compute Theta component
         this->transformST2Theta(rCurl.rTheta());

         // Compute Phi component
         this->transformST2Phi(rCurl.rPhi());
      }
   }

   template <typename TSpectralSHTraits, typename TSHTraits> void TorPolSpectralSHTransform<TSpectralSHTraits, TSHTraits>::initRTPDomains(SmartTruncation &rTrunc) const
   {
      // Fill in radial grid
      rTrunc->setRadialGrid(this->radBasis().at(0).pGrid());

      // Fill in theta grid (and cos, sin)
      rTrunc->setThetaGrid(this->hozBasis().at(0).pTheta(), this->hozBasis().at(0).pCosTheta(), this->hozBasis().at(0).pSinTheta());
   }

   template <typename TSpectralSHTraits, typename TSHTraits> void TorPolSpectralSHTransform<TSpectralSHTraits, TSHTraits>::transformPol2R(RTPScalar &rRComp)
   {
      // Extract real space values from the SH expansion
      SHTransform::transformSpec2RTP(rRComp);
   }

   template <typename TSpectralSHTraits, typename TSHTraits> void TorPolSpectralSHTransform<TSpectralSHTraits, TSHTraits>::transformPol2R(RTPScalar &rRComp, const typename TorPolSpectralSHTransform<TSpectralSHTraits, TSHTraits>::SpectralSHScalarType &polComp)
   {
      if(this->sshManipulator().atSendStage())
      {
         // Get temporary storage
         FDSHBackwardType &rLTmp = this->sshManipulator().provideBTmp();

         // Perform the backward Worland transform for Poloidal component: Q part
         this->spectralT().template multL<SetProduct>(rLTmp, polComp, &PolynomialType::projPol2Q);

         // Transpose the FDSHForwardType values to a suitable order for the SHTransform
         this->sshManipulator().send(rLTmp);
         this->sshManipulator().releaseTmp(rLTmp);
      } 

      if(this->sshManipulator().atRecvStage())
      {
         // Get temporary storage and receive data
         FDSHForwardType &rMTmp = this->sshManipulator().provideFTmp();
         this->sshManipulator().receive(rMTmp);

         // Extract real space values from the SH expansion
         SHTransform::transformSpec2RTP(rRComp, rMTmp);

         // free temporary storage
         this->sshManipulator().freeTmp(rMTmp);
      }
   }

}

#endif // TORPOLSPECTRALSHTRANSFORM_HPP
