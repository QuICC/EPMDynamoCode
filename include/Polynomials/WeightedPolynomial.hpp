/** \file WeightedPolynomial.hpp
 *  \brief Implementation of a general weighted polynomial
 */

#ifndef WEIGHTEDPOLYNOMIAL_HPP
#define WEIGHTEDPOLYNOMIAL_HPP

// System includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Polynomials/Polynomial.hpp"

namespace EPMDynamo {

   /**
    * @brief This class extends the Polynomial class to the weighted polynomials case
    *
    * This class extends the Polynomial class to the use of weighted polynomials. This
    * will be used as a component for the bases used in the Quadrature integrations.
    *
    * This is in principle exactly the same as the Polynomial class it simply adds the
    * quadrature weighted values of the polynomials to its attributes. The storage for the
    * weighted values is of the same kind than for the Polynomial class. 
    */
   class WeightedPolynomial : public Polynomial
   {
      public:

         /**
          * @brief Constructs the weighted polynomial from scratch
          *
          * @param param   Value of the fixed parameter
          * @param nX      Size of the grid
          * @param nN      Size of the varying parameter space
          * @param maxDiff Highest derivative order. It has a default value of 1
          */
         WeightedPolynomial(const int param, const int nX, const int nN, const int maxDiff = 1);

         /**
          * @brief Constructs the weighted polynomial on the given values
          *
          * The WeightedPolynomial is constructed on the given grid and quadrature weights
          *
          * @param param   Value of the fixed parameter
          * @param grid    SmartArray of the grid
          * @param nN      Size of the parameter space of the second parameter
          * @param weights SmartArray of the weights
          * @param maxDiff Highest derivative order. It has a default value of 1
          */
         WeightedPolynomial(const int param, const SmartArray grid, const int nN, const SmartArray weights, const int maxDiff = 1);

         /**
          * @brief Simple empty destructor
          */
         virtual ~WeightedPolynomial() {};

         /**
          * @brief Get polynomial weights
          */
         const Array& weights() const;

         /**
          * @brief Get SmartArry pointer to the weights
          */
         const SmartArray pWeights() const;

         /**
          * @brief Get weighted polynomial values
          */
         const Matrix& wPoly() const;

         /**
          * @brief Get weighted derivatives values
          */
         const Matrix& wDiff(const int i) const;

      protected:

         /**
          * @brief Set polynomial weights
          */
         Array& rWeights();

         /**
          * @brief Set weighted polynomial values
          */
         Matrix& rWPoly();

         /**
          * @brief Set weighted derivatives values
          *
          * @param i Ith derivative
          */
         Matrix& rWDiff(const int i);

         /**
          * @brief Initialise the weighted polynomials from scratch
          */
         virtual void initFull() = 0;

         /**
          * @brief Initialise the weighted polynomials with provided grid and weights
          */
         virtual void initPartial() = 0;

         /**
          * @brief Compute the weighted polynomials
          */
         virtual void computeWPolynomials() = 0;

         /**
          * @brief Compute the weighted polynomials
          */
         virtual void computeWDerivatives() = 0;

      private:

         /**
          * @brief SmartArray of the weights
          */
         SmartArray    mpWeights;

         /**
          * @brief Matrix of the weighted polynomials on all grid points
          */
         Matrix   mWPoly;

         /**
          * @brief Vector of Matrices of the weighted derivatives values on all grid points.
          */
         std::vector<Matrix>   mWDiff;

         /**
          * @brief Create vector of Matrices for the weighted derivatives
          */
         void initWDiff();

   };

   // Public inline functions
   // 

   inline const Array& WeightedPolynomial::weights() const
   {
      return *this->mpWeights;
   }

   inline const SmartArray WeightedPolynomial::pWeights() const
   {
      return this->mpWeights;
   }

   inline const Matrix& WeightedPolynomial::wPoly() const
   {
      return this->mWPoly;
   }

   inline const Matrix& WeightedPolynomial::wDiff(const int i) const
   {
      return this->mWDiff.at(i-1);
   }

   // Protected inline functions
   //

   inline Array& WeightedPolynomial::rWeights()
   {
      return *this->mpWeights;
   }

   inline Matrix& WeightedPolynomial::rWPoly()
   {
      return this->mWPoly;
   }

   inline Matrix& WeightedPolynomial::rWDiff(const int i)
   {
      return this->mWDiff.at(i-1);
   }

}

#endif // WEIGHTEDPOLYNOMIAL_HPP
