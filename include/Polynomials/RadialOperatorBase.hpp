/** \file RadialOperatorBase.hpp
 *  \brief Base of the implementation of a set of radial operators
 */

#ifndef RADIALOPERATORBASE_HPP
#define RADIALOPERATORBASE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "GridGenerators/LegendreGrid.hpp"
#include "Simulations/SimulationConstants.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the implementation of a set of radial operators
    *
    * \tparam TPolynomial Type of the underlying polynomial
    *
    * \epmTodo Review exactly how the energy integrals are computed
    */
   template <typename TPolynomial> class RadialOperatorBase: public TPolynomial
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
         RadialOperatorBase(const int pL, const int nX, const int nN);

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
         RadialOperatorBase(const int pL, const SmartArray grid, const int nN, const SmartArray weights);

         /**
          * @brief Simple empty destructor
          */
         virtual ~RadialOperatorBase() {};

         /**
          * @brief Get \f$l(l+1)\f$ factor
          */
         EPMFloat ll() const;

         /**
          * @brief Get \f$\frac{1}{l(l+1)}\f$ factor
          */
         EPMFloat ll_1() const;

         /**
          * @brief Get \f$\sqrt{l(l+1)}\f$ factor
          */
         EPMFloat sll() const;

         /**
          * @brief Get \f$\frac{1}{\sqrt{l(l+1)}}\f$ factor
          */
         EPMFloat sll_1() const;

         /**
          * @brief Get the energy matrix for a simple product
          */
         const Matrix&  eWeights() const;

      protected:

         /**
          * @brief Storage for different \f$l(l+1)\f$ factors
          *
          * first element: \f$l(l+1)\f$
          * second element: \f$\frac{1}{l(l+1)}\f$
          * third element: \f$\sqrt{l(l+1)}\f$
          * fourth element: \f$\frac{1}{\sqrt{l(l+1)}}\f$
          */
         Array    mLlFactors;

         /**
          * @brief Storage for \f$1/r\f$ values for efficiency
          */
         Array mR_1;

         /**
          * @brief Storage for the operator: \f$\frac{1}{r}\frac{d}{d r}(r P)\f$
          */
         Matrix   mR_1DRPoly;

         /**
          * @brief Store the laplacian projector operator
          */
         Matrix   mLaplacianProj;

         /**
          * @brief Store the laplacian spectral operator
          */
         Matrix   mLaplacianSpec;

         /**
          * @brief Store the Energy integral weights of the polynomials
          */
         Matrix    mEWeights;

      private:
         /**
          * @brief Compute the "l(l+)" factors
          */
         void initLlFactors();

         /**
          * @brief Compute the stored 1/r array
          */

         void computeR_1();

         /**
          * @brief Compute the \f$\frac{1}{r}\frac{d}{d r} r P\f$ operators
          *          for the standard polynomials and the QS polynomials
          */
         void computeR_1DRPOperators();

         /**
          * @brief Compute the projector and spectral laplacian operators
          */
         void computeLaplacians();

         /**
          * @brief Compute the energy integral factors
          *
          * This computes and stores the values \f$ \int_{0}^{1} r^2 P d r\f$
          */
         void computeEnergyIntegrals();

         /**
          * @brief Correct operators computations in special cases (for example CSCS ouput)
          */
         void correctOperators();
   };

   //
   // Inline definitions of access of Ll factors values
   //

   template <typename TPolynomial> inline EPMFloat  RadialOperatorBase<TPolynomial>::ll() const
   {
      return this->mLlFactors(0);
   }

   template <typename TPolynomial> inline EPMFloat  RadialOperatorBase<TPolynomial>::ll_1() const
   {
      return this->mLlFactors(1);
   }

   template <typename TPolynomial> inline EPMFloat  RadialOperatorBase<TPolynomial>::sll() const
   {
      return this->mLlFactors(2);
   }

   template <typename TPolynomial> inline EPMFloat  RadialOperatorBase<TPolynomial>::sll_1() const
   {
      return this->mLlFactors(3);
   }

   //
   // Inline definitions of the spectral operators
   //

   template <typename TPolynomial> inline const Matrix& RadialOperatorBase<TPolynomial>::eWeights() const
   {
      return this->mEWeights;
   }

   template <typename TPolynomial> RadialOperatorBase<TPolynomial>::RadialOperatorBase(const int pL, const int nX, const int nN)
      : TPolynomial(pL, nX, nN), mLlFactors(4), mR_1(nX), mR_1DRPoly(nN, nX), mLaplacianProj(nN, nX), mLaplacianSpec(nN, nN), mEWeights(nN, nN)
   {
      // Precompute "l(l+1)" factors
      this->initLlFactors();

      // Compute 1/r values
      this->computeR_1();

      // Compute the two R_1DR P operators
      this->computeR_1DRPOperators();

      // Compute the laplacians
      this->computeLaplacians();

      // Compute the energy integrals
      this->computeEnergyIntegrals();

      // Correct operators in special cases (for example CSCS output)
      this->correctOperators();
   }

   template <typename TPolynomial> RadialOperatorBase<TPolynomial>::RadialOperatorBase(const int pL, const SmartArray grid, const int nN, const SmartArray weights)
      : TPolynomial(pL, grid, nN, weights), mLlFactors(4), mR_1(this->gridN()), mR_1DRPoly(nN, this->gridN()), mLaplacianProj(nN, this->gridN()), mLaplacianSpec(nN, nN), mEWeights(nN, nN)
   {
      // Precompute "l(l+1)" factors
      this->initLlFactors();

      // Compute 1/r values
      this->computeR_1();

      // Compute the two R_1DR P operators
      this->computeR_1DRPOperators();

      // Compute the laplacians
      this->computeLaplacians();

      // Compute the energy integrals
      this->computeEnergyIntegrals();

      // Correct operators in special cases (for example CSCS output)
      this->correctOperators();
   }

   template <typename TPolynomial> void RadialOperatorBase<TPolynomial>::initLlFactors()
   {
      this->mLlFactors(0) = static_cast<EPMFloat>(this->l()*(this->l()+1));
      this->mLlFactors(2) = std::sqrt(this->ll());

      if(this->ll() == 0.0)
      {
         this->mLlFactors(1) = 0.0;
         this->mLlFactors(3) = 0.0;
      } else
      {
         this->mLlFactors(1) = 1.0/this->ll();
         this->mLlFactors(3) = 1.0/this->sll();
      }
   }

   template <typename TPolynomial> void RadialOperatorBase<TPolynomial>::computeR_1()
   {
      for(int i=0; i < this->gridN(); ++i)
      {
         this->mR_1(i) = 1.0/this->grid()(i);
      }
   }

   template <typename TPolynomial> void RadialOperatorBase<TPolynomial>::computeR_1DRPOperators()
   {
      // Compute the operator of an the standard polynomial
      this->mR_1DRPoly = this->diff(1) + this->poly()*this->mR_1.asDiagonal();
   }

   template <typename TPolynomial> void RadialOperatorBase<TPolynomial>::computeLaplacians()
   {
      Matrix   tmp1(this->polyN(), this->gridN());
      Array rfactor1 = static_cast<EPMFloat>(4*this->l()+6)*this->grid().array().pow(this->l());
      Array rfactor2 = 4.0*this->grid().array().pow(this->l()+2);

      this->computeDerivativeBase(this->mLaplacianProj, 1);

      this->computeDerivativeBase(tmp1, 2);

      for(int n = 1; n < this->polyN(); ++n)
      {
         this->mLaplacianProj.row(n).array() *= rfactor1.transpose().array();
         this->mLaplacianProj.row(n).array() += rfactor2.transpose().array()*tmp1.row(n).array();

         // normalise the laplacian projector
         this->mLaplacianProj.row(n).array() *= this->normalisation()(n);
      }

      this->mLaplacianSpec = (this->mLaplacianProj*this->wPoly()).transpose().template triangularView<Eigen::StrictlyUpper>();
   }

   template <typename TPolynomial> void RadialOperatorBase<TPolynomial>::computeEnergyIntegrals()
   {
      int ptsN = 2*this->polyN() + this->l() + 1;

      // Initialise the grid generator
      LegendreGrid  generator(ptsN);

      // Storage for the temporary grid
      SmartArray eGrid(new Array(ptsN));
      SmartArray eWeights(new Array(ptsN));

      // get computed grid
      (*eGrid) = generator.grid();

      // get computed weights
      (*eWeights) = generator.weights();

      // Adapt weights to stretched interval
      *eWeights = (*eWeights)/2.0;

      // Strech grid points to 0, 1 interval
      eGrid->array() += 1.0;
      *eGrid = (*eGrid)/2.0;

      TPolynomial  tmpPoly(this->l(), eGrid, this->polyN(), eWeights); 

      Matrix   tmp(this->polyN(), ptsN);
      EPMFloat correctNorm;

      for(int n=0; n < this->polyN(); ++n)
      {
         for(int i=0; i < this->polyN(); ++i)
         {
            correctNorm = this->normalisation()(i)*this->normalisation()(n)/(tmpPoly.normalisation()(i)*tmpPoly.normalisation()(n))
            this->mEWeights(i,n) = (tmpPoly.poly().row(n).array()*tmpPoly.poly().row(i).array()).matrix().dot ((eGrid->array().pow(2)*eWeights->array()).matrix())*correctNorm;
         }
      }
   }

   template <typename TPolynomial> void RadialOperatorBase<TPolynomial>::correctOperators()
   {
      // Treat the CSCS case differently (has 2 extraneous grid points)
      if(SimulationConstants::isCSCSGrid())
      {
         // Set 1/r to zero for r=0
         this->mR_1(0) = 0.0;

         // Set special values for r=0 for l = 0
         if(this->l() == 0)
         {
            this->mR_1DRPoly.col(0).setConstant(0.0);
         }

         // Set special values for r=0 for l = 1
         if(this->l() == 1)
         {
            this->mR_1DRPoly.col(0) = 2.0*this->diff(1).col(0);
         }

         // Set special values for r=0 for l > 1
         if(this->l() > 1)
         {
            this->mR_1DRPoly.col(0).setConstant(0.0);
         }
      }
   }

}

#endif // RADIALOPERATORBASE_HPP
