/** \file ChebyshevPolynomial.hpp
 *  \brief Implementation of the Chebyshev polynomials
 */

#ifndef CHEBYSHEVPOLYNOMIAL_HPP
#define CHEBYSHEVPOLYNOMIAL_HPP

// System includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Polynomials/WeightedPolynomial.hpp"

namespace EPMDynamo {

   /**
    * @brief This class defines the Chebyshev polynomials \f$T_n\f$
    *
    * The Polynomials are not normalised.
    *
    * \bug This code needs to be tested
    * \bug Correct documentation, it still holds the Jacobi documentation!
    */
   class ChebyshevPolynomial: public WeightedPolynomial
   {
      public:

         /**
          * @brief Constructs the polynomial from scratch
          *
          * This constructor computes the values of the grid and the weights and then
          * compute the polynomial values based on these
          *
          * @param pL      Polynomial L order
          * @param nX      Size of the grid
          * @param nN      Size of polynomial basis
          */
         ChebyshevPolynomial(const int pL, const int nX, const int nN);

         /**
          * @brief Constructs the polynomial based on the given grid points and weights
          *
          * The values of the polynomial are computed on the given grid and weighted with
          * the given weights.
          *
          * @param pL      Polynomial L order
          * @param grid    Provided grid array
          * @param nN      Size of polynomial basis
          * @param weights Provided weights array
          */
         ChebyshevPolynomial(const int pL, const SmartArray grid, const int nN, const SmartArray weights);

         /**
          * @brief Simple empty destructor
          */
         virtual ~ChebyshevPolynomial() {};

         /**
          * @brief Get the boundary values polynomials
          */
         const Array&   bpoly() const;

         /**
          * @brief Get the boundary values derivatives
          *
          * @param i Ith derivative
          */
         const Array&   bdiff(const int i) const;

      protected:
         /**
          * @brief Static value for the highest implemented derivative
          */
         static const int  msMaxDiff = 2;

         /**
          * @brief Complete initialisation of the polynomials
          */
         virtual void initFull();

         /**
          * @brief Initialise the polynomials from the given grid and weights
          */
         virtual void initPartial();

         /**
          * @brief Compute the grid points
          *
          * The grid points used are the Chebyshev points given by
          * \f$r_i = \cos\left(\frac{\pi}{4}\frac{(2i -1)}{N}\right)\f$ where \f$N\f$ is
          * the total number of grid points. The standard chebyshev grid is given by
          * \f$x_i = \cos\left(\frac{\pi}{2}\frac{(2i -1)}{N}\right)\f$ which is between
          * \f$-1\f$ and \f$1\f$. But we are interessed in a radial grid from \f$0\f$ to
          * \f$1\f$ which is obtain though \f$x = 2 r^2 -1 \f$.
          *
          * The used weights are the chebyshev weights given by \f$ w_i = \frac{\pi}{N}\f$
          * where \f$N\f$ is the total number of weights.
          */
         void computeWeightedGrid();

         /**
          * @brief Compute the Polynomials
          *
          * The polynomials are computed through their recurrence relation. 
          */
         virtual void computePolynomials();

         /**
          * @brief Compute the Polynomials on the boundary
          *
          * The polynomials are computed through their recurrence relation. 
          */
         virtual void computeBoundaryPolynomials();

         /**
          * @brief Compute the Polynomials
          */
         virtual void computeWPolynomials();

         /**
          * @brief Compute the Derivatives
          *
          * \warning This routine currently only computes the first derivative
          */
         virtual void computeDerivatives();

         /**
          * @brief Compute the Derivatives on the boundary
          *
          * \warning This routine currently only computes the first and second derivatives
          */
         virtual void computeBoundaryDerivatives();

         /**
          * @brief Compute the Derivatives
          */
         virtual void computeWDerivatives();

         /**
          * @brief Set the values of the boundary polynomials
          */
         Array& rBPoly();

         /**
          * @brief Set the boundary values derivatives
          *
          * @param i Ith derivative
          */
         Array&   rBDiff(const int i);

      private:
         /**
          * @brief Values of the polynomials on the boundary
          */
         Array mBPoly;

         /**
          * @brief Values of the derivatives of the polynomials on the boundary
          */
         std::vector<Array>  mBDiff;

         /**
          * @brief Initialise vector of boundary derivatives
          */
         void initBDiff();

         /**
          * @brief Compute the first derivative
          *
          * The derivative is also computed through a recurrence relation.
          * \f[
          *    \frac{d P_n^{(\alpha,\beta)}(2r^2 -1)}{d r} = 2r(n +\alpha+\beta+1) P_{n-1}^{(\alpha+1,\beta+1)}(2r^2-1)
          * \f]
          */
         void computeFirstDerivative();

         /**
          * @brief Compute the first derivative on the boundary
          *
          * The derivative is also computed through a recurrence relation.
          * \f[
          *    \frac{d P_n^{(\alpha,\beta)}(2r^2 -1)}{d r} = 2r(n +\alpha+\beta+1) P_{n-1}^{(\alpha+1,\beta+1)}(2r^2-1)
          * \f]
          */
         void computeBoundaryFirstDerivative();

         /**
          * @brief Compute the second derivative
          *
          * The derivative is also computed through a recurrence relation.
          * \f{eqnarray*}{
          *    \frac{d^2 P_n^{(\alpha,\beta)}(2r^2 -1)}{d r^2} & = & 4r^2(n +\alpha+\beta+1)(n +\alpha+\beta+2) P_{n-2}^{(\alpha+2,\beta+2)}(2r^2-1)\\
          *       & & + 2(n +\alpha+\beta+1) P_{n-1}^{(\alpha+1,\beta+1)}(2r^2-1)
          * \f}
          */
         void computeSecondDerivative();

         /**
          * @brief Compute the second derivative on the boundary
          *
          * The derivative is also computed through a recurrence relation.
          * \f{eqnarray*}{
          *    \frac{d^2 P_n^{(\alpha,\beta)}(2r^2 -1)}{d r^2} & = & 4r^2(n +\alpha+\beta+1)(n +\alpha+\beta+2) P_{n-2}^{(\alpha+2,\beta+2)}(2r^2-1)\\
          *       & & + 2(n +\alpha+\beta+1) P_{n-1}^{(\alpha+1,\beta+1)}(2r^2-1)
          * \f}
          */
         void computeBoundarySecondDerivative();

         /**
          * @brief Get recurrence coefficients \f$a_n = \frac{(2n+\alpha+\beta-1)(\alpha^2-\beta^2)}{2n(n+\alpha+\beta)(2n+\alpha+\beta-2)}\f$ 
          *
          * This routines use the stored \f$\alpha\f$ and \f$\beta\f$ values
          *
          * @param n Order n of the polynomial
          */
         EPMFloat recurrenceA(const int n) const;

         /**
          * @brief Get recurrence coefficients \f$b_n = \frac{(2n+\alpha+\beta-1)(2n +\alpha+\beta)}{2n(n+\alpha+\beta)}\f$ 
          *
          * This routines use the stored \f$\alpha\f$ and \f$\beta\f$ values
          *
          * @param n Order n of the polynomial
          */
         EPMFloat recurrenceB(const int n) const;

         /**
          * @brief Get recurrence coefficients \f$c_n = \frac{(n+\alpha-1)(n+\beta-1)(2n+\alpha+\beta)}{n(n+\alpha+\beta)(2n+\alpha+\beta-2)} \f$ 
          *
          * This routines use the stored \f$\alpha\f$ and \f$\beta\f$ values
          *
          * @param n Order n of the polynomial
          */
         EPMFloat recurrenceC(const int n) const;

         /**
          * @brief Compute normalisation factor required for weight
          *
          * The normalisation factor is computed with the help of the
          * quadrature rule.
          *
          * @param n Order of the polynomial
          */
         EPMFloat normaliseW(const int n);
   };

   // public inline functions
   //
   inline const Array& ChebyshevPolynomial::bpoly() const
   {
      return this->mBPoly;
   }

   inline Array& ChebyshevPolynomial::rBPoly()
   {
      return this->mBPoly;
   }

   inline const Array& ChebyshevPolynomial::bdiff(const int i) const
   {
      return this->mBDiff.at(i-1);
   }

   inline Array& ChebyshevPolynomial::rBDiff(const int i)
   {
      return this->mBDiff.at(i-1);
   }
}

#endif // CHEBYSHEVPOLYNOMIAL_HPP
