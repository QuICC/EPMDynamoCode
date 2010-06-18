/** \file TorPolRadialOperator.hpp
 *  \brief Implementation of a general set of Toroidal/Poloidal radial operators
 */

#ifndef TORPOLRADIALOPERATOR_HPP
#define TORPOLRADIALOPERATOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Polynomials/TorPolRadialOperatorBase.hpp"
#include "Polynomials/PolynomialOperator.hpp"
#include "Polynomials/PolyOps/PolyOpsTypes.hpp"

namespace EPMDynamo {

   /**
    * @brief This class defines the TorPolRadialOperator
    *
    * This is a extension to the polynomial given as template parameter
    * which allows to use it as a radial basis. It defines and creates
    * the required special operators.
    *
    * \tparam TPolynomial Type of the underlying polynomial
    *
    * \bug Rewrite the polynomial operators and simpler and cleaner way
    */
   template <typename TPolynomial> class TorPolRadialOperator: public TorPolRadialOperatorBase<TPolynomial>
   {

      public:
         /**
          * @brief Constructs the polynomials and sets up addition operators
          *
          * This constructor computes the values of the grid and the weights and then
          * compute the polynomial values based on these
          *
          * @param pL      Polynomial L order
          * @param nX      Size of the grid
          * @param nN      Size of polynomial basis
          */
         TorPolRadialOperator(const int pL, const int nX, const int nN);

         /**
          * @brief Constructs the polynomials with given points and weighs and sets up
          * addition operators
          *
          * The values of the polynomial are computed on the given grid and weighted with
          * the given weights.
          *
          * @param pL      Polynomial L order
          * @param grid    Provided grid array
          * @param nN      Size of polynomial basis
          * @param weights Provided weights array
          */
         TorPolRadialOperator(const int pL, const SmartArray grid, const int nN, const SmartArray weights);

         /**
          * @brief Need copy constructor to avoid shallow copy!
          *
          * @param poly Polynomial to copy
          */
         TorPolRadialOperator(const TorPolRadialOperator &poly);

         /**
          * @brief Need copy assignement to avoid shallow copy!
          *
          * @param poly Polynomial to assign
          */
         TorPolRadialOperator&  operator=(const TorPolRadialOperator &poly);

         /**
          * @brief Simple empty destructor
          */
         virtual ~TorPolRadialOperator() {};

         //
         // Integration operators
         //
         
         /**
          * @brief Get the integration operator
          *
          * This is the integration operator for a simple scalar
          */
         const PolynomialOperator<POperator>&  intg() const;

         /**
          * @brief Get the integration operator from a QST T component to the Toroidal
          *          component
          *
          * \todo Add more information
          */
         const PolynomialOperator<CPOperator<DynamoFloat> >&  intgT2Tor() const;

         /**
          * @brief Get the integration operator from a QST Q component to the Poloidal
          *          component
          *
          * \todo Add more information
          */
         const PolynomialOperator<CLPOperator<DynamoFloat> >&  intgQ2Pol() const;
         
         /**
          * @brief Get the integration operator for the Curl projection
          *
          * \todo Add more information
          */
         const PolynomialOperator<CPOperator<DynamoFloat> >&  intgT2CurlProj() const;
         
         /**
          * @brief Get the integration operator for the CurlCurl projection
          *
          * \todo Add more information
          */
         const PolynomialOperator<LPOperator>&  intgQ2CurlCurlProj() const;
         
         /**
          * @brief Get the integration operator for the CurlCurl projection
          *
          * \todo Add more information
          */
         const PolynomialOperator<POperator>&  intgS2CurlCurlProj() const;

         //
         // Projection operators
         //

         /**
          * @brief Get the projection operator for a standard \f$W_n^{l}\f$ expansion
          *
          * \todo Add more information
          */
         const PolynomialOperator<POperator>&  proj() const;

         /**
          * @brief Get the projection operator of the first derivative of a 
          *          \f$W_n^{l-1}\f$ expansion
          *
          * \todo Add more information
          */
         const PolynomialOperator<POperator>&  dProj() const;

         /**
          * @brief Get the projection operator of the second derivative of a 
          *          \f$W_n^{l-1}\f$ expansion
          *
          * \todo Add more information
          */
         const PolynomialOperator<POperator>&  d2Proj() const;

         /**
          * @brief Get the projection operator from a spectral scalar to the theta
          *          or phi component of the gradient
          *
          * \todo Add more information
          */
         const PolynomialOperator<RPOperator>&  proj2GradTP() const;

         /**
          * @brief Get the projection operator from a Toroidal component to the Q
          * component of the curl
          *
          * \todo Add more information
          */
         const PolynomialOperator<CRPOperator<DynamoFloat> >&  projTor2CurlQ() const;

         /**
          * @brief Get the projection operator from a Toroidal component to the S
          * component of the curl
          *
          * \todo Add more information
          */
         const PolynomialOperator<CPOperator<DynamoFloat> >&  projTor2CurlS() const;

         /**
          * @brief Get the projection operator from a Poloidal component to the T
          * component of the curl
          *
          * \todo Add more information
          */
         const PolynomialOperator<CPOperator<DynamoFloat> >&  projPol2CurlT() const;

         /**
          * @brief Get the projection operator from the Poloidal component to the QST Q 
          *          component
          *
          * \todo Add more information
          */
         const PolynomialOperator<CRPOperator<DynamoFloat> >&  projPol2Q() const;

         /**
          * @brief Get the projection operator from the Poloidal component to the QST S 
          *          component
          *
          * \todo Add more information
          */
         const PolynomialOperator<CPOperator<DynamoFloat> >&  projPol2S() const;

         /**
          * @brief Get the projection operator from the Toroidal component to the QST T 
          *          component
          *
          * \todo Add more information
          */
         const PolynomialOperator<CPOperator<DynamoFloat> >&  projTor2T() const;

         //
         // Spectral space operators
         //

         /**
          * @brief Get the spectral operator for the laplacian 
          *
          * \todo Add more information
          */
         const PolynomialOperator<POperator>&  specLaplacian() const;

      protected:

      private:
         /**
          * @brief Stored smart P operator
          */
         SmartPOp    mpIntg;

         /**
          * @brief Stored smart DP operator
          */
         SmartDPOp   mpIntgT2Tor;

         /**
          * @brief Stored smart DLP operator
          */
         SmartDLPOp  mpIntgQ2Pol;

         /**
          * @brief Stored smart DP operator
          */
         SmartDPOp   mpIntgT2CurlProj;

         /**
          * @brief Stored smart LP operator
          */
         SmartLPOp   mpIntgQ2CurlCurlProj;

         /**
          * @brief Stored smart P operator
          */
         SmartPOp    mpIntgS2CurlCurlProj;

         /**
          * @brief Stored smart P operator
          */
         SmartPOp    mpProj;

         /**
          * @brief Stored smart P operator
          */
         SmartPOp    mpDProj;

         /**
          * @brief Stored smart P operator
          */
         SmartPOp    mpD2Proj;

         /**
          * @brief Stored smart RP operator
          */
         SmartRPOp   mpProj2GradTP;

         /**
          * @brief Stored smart DRP operator
          */
         SmartDRPOp  mpProjTor2CurlQ;

         /**
          * @brief Stored smart DP operator
          */
         SmartDPOp   mpProjTor2CurlS;

         /**
          * @brief Stored smart DP operator
          */
         SmartDPOp   mpProjPol2CurlT;

         /**
          * @brief Stored smart DRP operator
          */
         SmartDRPOp  mpProjPol2Q;

         /**
          * @brief Stored smart DP operator
          */
         SmartDPOp   mpProjPol2S;

         /**
          * @brief Stored smart DP operator
          */
         SmartDPOp   mpProjTor2T;

         /**
          * @brief Stored smart P operator
          */
         SmartPOp    mpSpecLaplacian;

         /**
          * @brief Init all the PolynomialOperator<*> objects
          */
         void initPOperators();
   };

   //
   // Inline definitions of the integrators
   //

   template <typename TPolynomial> inline const PolynomialOperator<POperator>& TorPolRadialOperator<TPolynomial>::intg() const
   {
      return *this->mpIntg;
   }

   template <typename TPolynomial> inline const PolynomialOperator<CPOperator<DynamoFloat> >& TorPolRadialOperator<TPolynomial>::intgT2Tor() const
   {
      return *this->mpIntgT2Tor;
   }

   template <typename TPolynomial> inline const PolynomialOperator<CLPOperator<DynamoFloat> >& TorPolRadialOperator<TPolynomial>::intgQ2Pol() const
   {
      return *this->mpIntgQ2Pol;
   }

   template <typename TPolynomial> inline const PolynomialOperator<CPOperator<DynamoFloat> >& TorPolRadialOperator<TPolynomial>::intgT2CurlProj() const
   {
      return *this->mpIntgT2CurlProj;
   }

   template <typename TPolynomial> inline const PolynomialOperator<LPOperator>& TorPolRadialOperator<TPolynomial>::intgQ2CurlCurlProj() const
   {
      return *this->mpIntgQ2CurlCurlProj;
   }

   template <typename TPolynomial> inline const PolynomialOperator<POperator>& TorPolRadialOperator<TPolynomial>::intgS2CurlCurlProj() const
   {
      return *this->mpIntgS2CurlCurlProj;
   }

   //
   // Inline definitions of the projectors
   //

   template <typename TPolynomial> inline const PolynomialOperator<POperator>& TorPolRadialOperator<TPolynomial>::proj() const
   {
      return *this->mpProj;
   }

   template <typename TPolynomial> inline const PolynomialOperator<POperator>& TorPolRadialOperator<TPolynomial>::dProj() const
   {
      return *this->mpDProj;
   }

   template <typename TPolynomial> inline const PolynomialOperator<POperator>& TorPolRadialOperator<TPolynomial>::d2Proj() const
   {
      return *this->mpD2Proj;
   }

   template <typename TPolynomial> inline const PolynomialOperator<RPOperator>& TorPolRadialOperator<TPolynomial>::proj2GradTP() const
   {
      return *this->mpProj2GradTP;
   }

   template <typename TPolynomial> inline const PolynomialOperator<CRPOperator<DynamoFloat> >& TorPolRadialOperator<TPolynomial>::projTor2CurlQ() const
   {
      return *this->mpProjTor2CurlQ;
   }

   template <typename TPolynomial> inline const PolynomialOperator<CPOperator<DynamoFloat> >& TorPolRadialOperator<TPolynomial>::projTor2CurlS() const
   {
      return *this->mpProjTor2CurlS;
   }

   template <typename TPolynomial> inline const PolynomialOperator<CPOperator<DynamoFloat> >& TorPolRadialOperator<TPolynomial>::projPol2CurlT() const
   {
      return *this->mpProjPol2CurlT;
   }

   template <typename TPolynomial> inline const PolynomialOperator<CRPOperator<DynamoFloat> >& TorPolRadialOperator<TPolynomial>::projPol2Q() const
   {
      return *this->mpProjPol2Q;
   }

   template <typename TPolynomial> inline const PolynomialOperator<CPOperator<DynamoFloat> >& TorPolRadialOperator<TPolynomial>::projPol2S() const
   {
      return *this->mpProjPol2S;
   }

   template <typename TPolynomial> inline const PolynomialOperator<CPOperator<DynamoFloat> >& TorPolRadialOperator<TPolynomial>::projTor2T() const
   {
      return *this->mpProjTor2T;
   }

   //
   // Inline definitions of the spectral operators
   //

   template <typename TPolynomial> inline const PolynomialOperator<POperator>& TorPolRadialOperator<TPolynomial>::specLaplacian() const
   {
      return *this->mpSpecLaplacian;
   }

   template <typename TPolynomial> TorPolRadialOperator<TPolynomial>::TorPolRadialOperator(const int pL, const int nX, const int nN)
      : TorPolRadialOperatorBase<TPolynomial>(pL, nX, nN)
   {
      // Initialise Polynomial operators
      this->initPOperators();
   }

   template <typename TPolynomial> TorPolRadialOperator<TPolynomial>::TorPolRadialOperator(const int pL, const SmartArray grid, const int nN, const SmartArray weights)
      : TorPolRadialOperatorBase<TPolynomial>(pL, grid, nN, weights)
   {
      // Initialise Polynomial operators
      this->initPOperators();
   }

   template <typename TPolynomial> TorPolRadialOperator<TPolynomial>::TorPolRadialOperator(const TorPolRadialOperator<TPolynomial> &other)
      : TorPolRadialOperatorBase<TPolynomial>(other)
   {
      // Initialise Polynomial operators
      this->initPOperators();
   }

   template <typename TPolynomial> TorPolRadialOperator<TPolynomial> & TorPolRadialOperator<TPolynomial>::operator=(const TorPolRadialOperator &other)
   {
      if(this == &other)
      {
         return *this;
      } else
      {
         // Copy object
         *this = other;

         // Recreate the smart pointers
         this->initPOperators();

         return *this;
      }
   }

   template <typename TPolynomial> void TorPolRadialOperator<TPolynomial>::initPOperators()
   {
      //
      this->mpIntg = SmartPOp(new PolynomialOperator<POperator>(this->wPoly()));

      //
      this->mpIntgT2Tor = SmartDPOp(new PolynomialOperator<CPOperator<DynamoFloat> >(-this->sll_1(), this->wPoly()));

      //
      this->mpIntgQ2Pol = SmartDLPOp(new PolynomialOperator<CLPOperator<DynamoFloat> >(this->ll_1(), this->grid(), this->wPoly()));

      //
      this->mpIntgT2CurlProj = SmartDPOp(new PolynomialOperator<CPOperator<DynamoFloat> >(-this->sll_1(), this->wPoly()));

      //
      this->mpIntgQ2CurlCurlProj = SmartLPOp(new PolynomialOperator<LPOperator>(this->mR_1, this->wPoly()));

      //
      this->mpIntgS2CurlCurlProj = SmartPOp(new PolynomialOperator<POperator>(this->mS2CurlCurlIntegrator));

      //
      this->mpProj = SmartPOp(new PolynomialOperator<POperator>(this->poly()));

      //
      this->mpDProj = SmartPOp(new PolynomialOperator<POperator>(this->diff(1)));

      //
      this->mpD2Proj = SmartPOp(new PolynomialOperator<POperator>(this->diff(2)));

      //
      this->mpProj2GradTP = SmartRPOp(new PolynomialOperator<RPOperator>(this->mR_1, this->poly()));

      //
      this->mpProjTor2CurlQ = SmartDRPOp(new PolynomialOperator<CRPOperator<DynamoFloat> >(this->ll(), this->mR_1, this->poly()));

      //
      this->mpProjTor2CurlS = SmartDPOp(new PolynomialOperator<CPOperator<DynamoFloat> >(this->sll(), this->mR_1DRPoly));

      //
      this->mpProjPol2CurlT = SmartDPOp(new PolynomialOperator<CPOperator<DynamoFloat> >(this->sll(), this->mLaplacianProj));

      //
      this->mpProjPol2Q = SmartDRPOp(new PolynomialOperator<CRPOperator<DynamoFloat> >(this->ll(), this->mR_1, this->poly()));

      //
      this->mpProjPol2S = SmartDPOp(new PolynomialOperator<CPOperator<DynamoFloat> >(this->sll(), this->mR_1DRPoly));

      //
      this->mpProjTor2T = SmartDPOp(new PolynomialOperator<CPOperator<DynamoFloat> >(-this->sll(), this->poly()));

      //
      this->mpSpecLaplacian = SmartPOp(new PolynomialOperator<POperator>(this->mLaplacianSpec));
   }

}

#endif // TORPOLRADIALOPERATOR_HPP
