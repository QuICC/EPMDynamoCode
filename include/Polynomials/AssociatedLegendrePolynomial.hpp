/** \file AssociatedLegendrePolynomial.hpp
 *  \brief Implementation of the Associated Legendre polynomial
 */

#ifndef ASSOCIATEDLEGENDREPOLYNOMIAL_HPP
#define ASSOCIATEDLEGENDREPOLYNOMIAL_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Polynomials/WeightedPolynomial.hpp"
#include "General/MathConstants.hpp"

namespace EPMDynamo {

   /**
    * @brief This class defines the Associated Legendre polynomials
    *
    * This class represents the Associated Legendre polynomials. They are computed for
    * a fixed values of \f$m\f$ of the polynomial \f$P_{l}^{m}(\cos\theta)\f$. This means
    * that the values stored in a AssociatedLegendrePolynomial object are all the \f$l\f$ modes for
    * a fixed \f$m\f$.
    *
    * The computed polynomials are Schmidt quasi-normalised.
    *
    * \epmBug Explain the factor 2 required for m=0.
    *
    * \epmBug Explain the factor \f$\sqrt{2}\f$ required.
    */
   class AssociatedLegendrePolynomial: public WeightedPolynomial
   {
      public:

         /**
          * @brief Constructs the polynomial from scratch.
          *
          * This constructors computes the values of the grid and the weights and then
          * computes the polynomial values based on these
          *
          * @param nX      Size of the grid
          * @param nL      Number of l modes
          * @param pM      Polynomial m degree
          */
         AssociatedLegendrePolynomial(const int pM, const int nX, const int nL);

         /**
          * @brief Constructs the polynomial based on the given grid points and weights
          *
          * The values of the polynomial are computed on the given grid and weighted with
          * the given weights.
          *
          * @param grid    Provided grid array
          * @param nL      Number of l modes
          * @param weights Provided weights array
          * @param pM      Polynomial m degree
          */
         AssociatedLegendrePolynomial(const int pM, const SmartArray grid, const int nL, const SmartArray weights);

         /**
          * @brief Simple empty destructor
          */
         virtual ~AssociatedLegendrePolynomial() {};

         /**
          * @brief  Get the harmonic degree of the polynomial
          */
         int   m() const;

         /**
          * @brief Get theta coordinates
          */
         const Array&    theta() const;

         /**
          * @brief Get Sin(theta) coordinates
          */
         const Array&    sinTheta() const;

         /**
          * @brief Get Cos(theta) coordinates
          */
         const Array&    cosTheta() const;

         /**
          * @brief Get pointer to theta coordinates
          */
         const SmartArray&    pTheta() const;

         /**
          * @brief Get pointer to Sin(theta) coordinates
          */
         const SmartArray&    pSinTheta() const;

         /**
          * @brief Get pointer to Cos(theta) coordinates
          */
         const SmartArray&    pCosTheta() const;

      protected:

         /**
          * @brief Complete initialisation of the polynomials
          */
         virtual void initFull();

         /**
          * @brief Initialise the polynomials from given grid and weights
          */
         virtual void initPartial();

         /**
          * @brief Compute the grid points and the corresponding weights
          *
          * The grid points and the weights are optained through solving a tridiagonal matrix constructed
          * from the three term recurrence relation (see Golub's paper) 
          */
         void computeWeightedGrid();

         /**
          * @brief Compute the theta grid points and other useful grid values
          */
         void computeThetaGrids();

         /**
          * @brief Compute the Polynomials
          */
         virtual void computePolynomials();

         /**
          * @brief Compute the weighted Polynomials
          */
         virtual void computeWPolynomials();

         /**
          * @brief Compute the Derivatives
          */
         virtual void computeDerivatives();

         /**
          * @brief Compute the Derivatives
          */
         virtual void computeWDerivatives();

         /**
          * @brief Set theta coordinates
          */
         Array&    rTheta();

         /**
          * @brief Set Sin(theta) coordinates
          */
         Array&    rSinTheta();

         /**
          * @brief Set Cos(theta) coordinates
          */
         Array&    rCosTheta();

      private:

         /**
          * @brief Static variable setting the highest implemented derivative
          */
         static const int  msMaxDiff = 1;

         /**
          * @brief Array holding the theta values
          */
         SmartArray    mpTheta;

         /**
          * @brief Array holding the sin(theta) values
          */
         SmartArray    mpSinTheta;

         /**
          * @brief Compute normalisation factor for the weights 
          *
          * The weights have also to include the normalisation correction due to the fact
          * that the polynomials are not orthonormal but Schmidt quasi normalised. This
          * is required to avoid modifiying the values during the transforms.
          *
          * @param l Order n of the polynomial
          */
         EPMFloat normaliseW(const int l) const;

         /**
          * @brief Get recurrence coefficient \f$a_n = 0\f$ 
          *
          * \warning This is not the recurrence relation for the associated legendre 
          * polynomial but for the "classical" legendre polynomial. This is required to
          * compute the grid points and the related weights for the quadrature rule.
          *
          * @param n Order n of the polynomial
          */
         EPMFloat recurrenceA(const int n) const;

         /**
          * @brief Get recurrence coefficient \f$b_n = \frac{2n -1}{n}\f$ 
          *
          * \warning This is not the recurrence relation for the associated legendre 
          * polynomial but for the "classical" legendre polynomial. This is required to
          * compute the grid points and the related weights for the quadrature rule.
          *
          * @param n Order n of the polynomial
          */
         EPMFloat recurrenceB(const int n) const;

         /**
          * @brief Get recurrence coefficient \f$c_n = \frac{n-1}{n}\f$ 
          *
          * \warning This is not the recurrence relation for the associated legendre 
          * polynomial but for the "classical" legendre polynomial. This is required to
          * compute the grid points and the related weights for the quadrature rule.
          *
          * @param n Order n of the polynomial
          */
         EPMFloat recurrenceC(const int n) const;

         /**
          * @brief Compute the theta coordinates
          */
         void computeTheta();

         /**
          * @brief Compute the sin(theta)
          */
         void computeSinTheta();

         /**
          * @brief Compute the Associate Legendre polynomial \f$P_m^m(\cos\theta)\f$
          *
          * This method computes the polynomial \f$P_m^m(\cos\theta)\f$ using as recurrence relation.
          * The recurrencre relation directly computes the Schmidt quasi-normalised value
          *
          * \f[
          *    P_m^m(\cos\theta) = (-1)^m \sqrt{\frac{(2m-1)!!}{(2m)!!}} (\sin(\theta))^m
          * \f]
          *
          * @param m Harmonic order m
          * @param idx Index
          */
         void computePmm(const int m, const int idx = 0);

         /**
          * @brief Compute the Associate Legendre polynomial \f$P_{m+1}^m(\cos\theta)\f$
          *
          * This method computes the polynomial \f$P_{m+1}^m(\cos\theta)\f$ using as recurrence relation.
          * The recurrencre relation directly computes the Schmidt quasi-normalised value
          *
          * \f[
          *    P_{m+1}^m(\cos\theta) = \sqrt{2m +1}\cos(\theta)P_m^m(\cos\theta)
          * \f]
          *
          * @param m Harmonic order m
          * @param idx Index
          */
         void computePmm1(const int m, const int idx = 1);

         /**
          * @brief Compute the Associate Legendre polynomial \f$P_l^m(\cos\theta)\f$
          *
          * This method computes the polynomial \f$P_{l}^m(\cos\theta)\f$ using as recurrence relation.
          * The recurrencre relation directly computes the Schmidt quasi-normalised value
          *
          * \f[
          *    P_l^m(\cos\theta) = \frac{1}{\sqrt{(l+m)(l-m)}}\left( (2l-1)\cos\theta P_{l-1}^m - \sqrt{(l-m-1)(l+m-1)} P_{l-2}^m \right)
          * \f]
          *
          * @param l Harmonic degree l
          * @param m Harmonic order m
          * @param idx Index
          */
         void computePlm(int l, const int m, int idx = 2);

         /**
          * @brief Compute the Associate Legendre polynomial \f$\frac{d}{d\theta}P_m^m(\cos\theta)\f$
          *
          * This method computes the polynomial \f$\frac{d}{d\theta}P_{m}^m(\cos\theta)\f$ using as recurrence relation.
          * The recurrencre relation directly computes the Schmidt quasi-normalised value
          *
          * \f[
          *    \frac{d}{d\theta}P_m^m(\cos\theta) = (-1)^m \sqrt{\frac{(2m-1)!!}{(2m)!!}} m \cos\theta(\sin(\theta))^{m-1}
          * \f]
          *
          * @param m Harmonic order m
          * @param idx Index
          */
         void computeDPmm(const int m, const int idx = 0);

         /**
          * @brief Compute the Associate Legendre polynomial \f$\frac{d}{d\theta}P_{m+1}^m(\cos\theta)\f$
          *
          * This method computes the polynomial \f$\frac{d}{d\theta}P_{m+1}^m(\cos\theta)\f$ using as recurrence relation.
          * The recurrencre relation directly computes the Schmidt quasi-normalised value
          *
          * \f[
          *    \frac{d}{d\theta}P_{m+1}^m(\cos\theta) = \sqrt{2m + 1}\left( \cos\theta \frac{d}{d\theta}P_{m}^m(\cos\theta) - \sin\theta P_m^m(\cos\theta)\right)
          * \f]
          *
          * @param m Harmonic order m
          * @param idx Index
          */
         void computeDPmm1(const int m, const int idx = 1);

         /**
          * @brief Compute the Associate Legendre polynomial \f$\frac{d}{d\theta}P_l^m(\cos\theta)\f$
          *
          * This method computes the polynomial \f$\frac{d}{d\theta}P_{l}^m(\cos\theta)\f$ using as recurrence relation.
          * The recurrencre relation directly computes the Schmidt quasi-normalised value
          *
          * \f[
          *    \frac{d}{d\theta}P_l^m(\cos\theta) = \frac{1}{\sqrt{(l+m)(l-m)}}\left( (2l -1)\left(\cos\theta \frac{d}{d\theta}P_{l-1}^m(\cos\theta)- \sin\theta P_{l-1}^{m}\right) - \sqrt{(l-m-1)(l+m-1)}\frac{d}{d\theta}P_{l-2}^m(\cos\theta) \right)
          * \f]
          *
          * @param l Harmonic degree l
          * @param m Harmonic order m
          * @param idx Index
          */
         void computeDPlm(int l, const int m, int idx = 2);

         /**
          * @brief Compute the polynomial basis
          */
         void computePlmBasis();

         /**
          * @brief Compute the polynomial derivative basis
          */
         void computeDPlmBasis();
   };

   // Public inline functions
   //
   inline int AssociatedLegendrePolynomial::m() const
   {
      return this->mParam;
   }

   inline const Array& AssociatedLegendrePolynomial::theta() const
   {
      return *this->mpTheta;
   }

   inline const Array& AssociatedLegendrePolynomial::sinTheta() const
   {
      return *this->mpSinTheta;
   }

   inline const Array& AssociatedLegendrePolynomial::cosTheta() const
   {
      return this->grid();
   }

   inline const SmartArray& AssociatedLegendrePolynomial::pTheta() const
   {
      return this->mpTheta;
   }

   inline const SmartArray& AssociatedLegendrePolynomial::pSinTheta() const
   {
      return this->mpSinTheta;
   }

   inline const SmartArray& AssociatedLegendrePolynomial::pCosTheta() const
   {
      return this->pGrid();
   }

   inline Array& AssociatedLegendrePolynomial::rTheta()
   {
      return *this->mpTheta;
   }

   inline Array& AssociatedLegendrePolynomial::rSinTheta()
   {
      return *this->mpSinTheta;
   }

   inline Array& AssociatedLegendrePolynomial::rCosTheta()
   {
      return this->rGrid();
   }

}

#endif // ASSOCIATEDLEGENDREPOLYNOMIAL_HPP
