/** \file WorlandPolynomial.hpp
 *  \brief Implementation of the Worland polynomial
 */

#ifndef WORLANDPOLYNOMIAL_HPP
#define WORLANDPOLYNOMIAL_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Polynomials/OnesidedJacobi.hpp"

namespace EPMDynamo {

   /**
    * @brief This class defines the WorlandPolynomial polynomials
    *
    * This class represents the WorlandPolynomial polynomials. They are computed for a fixed value
    * of \f$l\f$ of the polynomial
    * \f$W_n^l(2 r^2 -1) = r^l P_n^{(\alpha,\beta)}(2r^2-1)\f$, where \f$\alpha = -1/2\f$
    * and \f$\beta = l -1/2\f$. This means that the values stored in a WorlandPolynomial object are
    * all the \f$n\f$ orders for a fixed \f$l\f$.
    *
    * The polynomials are not normalised but the weighting is normalised, i.e. forward + backward transforms doesn't add an scaling.
    */
   class WorlandPolynomial: public OnesidedJacobi
   {

      public:

         /**
          * @brief Constructs the polynomials from scratch
          *
          * This constructor computes the values of the grid and the weights and then
          * compute the polynomial values based on these
          *
          * @param pL      Polynomial L order
          * @param nX      Size of the grid
          * @param nN      Size of polynomial basis
          */
         WorlandPolynomial(const int pL, const int nX, const int nN);

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
         WorlandPolynomial(const int pL, const SmartArray grid, const int nN, const SmartArray weights);

         /**
          * @brief Simple empty destructor
          */
         virtual ~WorlandPolynomial() {};

         /**
          * @brief Get the normalisation factors
          */
         const Array& normalisation() const; 

      protected:
         /**
          * @brief Complete initialisation of the polynomials
          */
         virtual void initFull();

         /**
          * @brief Initialise the polynomials from the given grid and weights
          */
         virtual void initPartial();

         /**
          * @brief Compute the Polynomials
          *
          * Once the OnesidedJacobi are built, computing the polynomials is a simple 
          * question of multipling by the respective \f$r^l\f$ factors.
          */
         virtual void computePolynomials();

         /**
          * @brief Compute the Polynomials
          */
         virtual void computeWPolynomials();

         /**
          * @brief Compute the Derivatives
          */
         virtual void computeDerivatives();

         /**
          * @brief Compute the Derivatives on the boundary
          */
         virtual void computeBoundaryDerivatives();

         /**
          * @brief Compute the Derivatives
          */
         virtual void computeWDerivatives();

         /**
          * @brief Storage for the normalisation factors
          */
         Array mNormalisation;

      private:
         /**
          * @brief Compute the first Derivatives
          */
         void computeFirstDerivative();

         /**
          * @brief Compute the first Derivatives on the  boundary
          */
         void computeBoundaryFirstDerivative();

         /**
          * @brief Compute the second Derivatives
          */
         void computeSecondDerivative();

         /**
          * @brief Compute the second Derivatives on the boundary
          */
         void computeBoundarySecondDerivative();

         /**
          * @brief Normalise all the polynomials
          */
         void normalise();

         /**
          * @brief Compute normalisation factor required for weight
          *
          * The normalisation has no analytical form unfortunately. The normalisation 
          * factor is computed with the help of the quadrature rule.
          *
          * @param n Order of the polynomial
          */
         EPMFloat inverseNorm(const int n);

         /**
          * @brief Correct the polynomials for special cases (for example CSCS output)
          */
         void correctPolynomials();
   };
}

#endif // WORLANDPOLYNOMIAL_HPP
