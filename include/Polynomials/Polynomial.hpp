/** \file Polynomial.hpp
 *  \brief Implementation of a general polynomial
 */

#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /** 
    * @brief A class representing a two parameters polynomial.
    *
    * The Polynomial class is an abstract representation of a two parameters polynomials 
    * on a discretised grid. To fix the notation and terms let's define such a poynomial 
    * P in the following way: \f$P_a^b(x) \f$. \f$a\f$ and \f$b\f$ are the two parameters
    * of the polynomial.
    *
    * The mathematical operations that are easily performed concern array and matrix
    * operations. But the values of the polynomial is a 3 dimensional array. To represent
    * it in memory, one of the parameters will be fixed, while all the possible
    * combination of the grid values and the second parameter will be stored in a Matrix.
    *
    * Which parameter will be fixed and which one will be varied over its value space
    * will depend on the polynomial and its use.
    *
    * The polynomial should be an orthogonal polynomial that can be represented with a
    * three term recurrence relation:
    *
    * \f$ P_n(x) = (a_n + b_n x)P_{n-1} - c_n P_{n-2}(x)\f$
    *
    * The polynomials values and the derivatives are stored with each row of the
    * matrix representing one value of the varying polynomial on the whole grid. This 
    * means that the fast direction through the array is by going through all values
    * of the paramter for each grid point. 
    *
    */
   class Polynomial
   {
      public:
         
         /** 
          * @brief Constructs the polynomial from scratch
          *
          * This constructor creates a polynomial from scratch. In a first step the 
          * discretised grid gets computed and then the values of the polynomial are
          * computed.
          *
          * @param param   Value of the fixed parameter of the polynomial
          * @param nX      Size of the grid
          * @param nN      Size of the parameter space of the second parameter
          * @param maxDiff Highest derivative order. It has a default value of 1
          */
         Polynomial(const int param, const int nX, const int nN, const int maxDiff = 1);

         /** 
          * @brief Constructs the polynomials on the given grid values
          *
          * This constructor builds the polynomial from the given grid values. Only a
          * pointer (Smart version) is stored to the grid.
          *
          * @param param   Value of the fixed parameter of the polynomial
          * @param grid    SmartArray pointer to the grid values
          * @param nN      Size of the parameter space of the varying parameter
          * @param maxDiff Highest derivative order. It has a default value of 1
          */
         Polynomial(const int param, const SmartArray grid, const int nN, const int maxDiff = 1);

         /**
          * @brief Empty destructor
          */
         virtual ~Polynomial() {};

         // Accessors
         //
         /**
          * @brief Get grid size
          */
         int gridN() const;

         /**
          * @brief Get the size of parameter space of the varying parameter
          */
         int polyN() const;

         /**
          * @brief Get grid points
          *
          * This routines returns a columnar array
          */
         const Array& grid() const;

         /**
          * @brief Get the smart pointer to the grid
          */
         const SmartArray& pGrid() const;

         /**
          * @brief Get polynomial values
          */
         const Matrix& poly() const;

         /**
          * @brief Get the i-th derivative values
          *
          * @param i    Order of the derivative
          */
         const Matrix& diff(const int i) const;

      protected:

         /**
          * @brief Value of the fixed parameter.
          */
         int mParam;

         /**
          * @brief Get highest derivatives order
          */
         int  maxDiff();

         /**
          * @brief Get the recurrence coefficient \f$a_n\f$
          *
          * @param n get Nth coefficient
          */
         virtual DynamoFloat recurrenceA(const int n) const = 0;

         /**
          * @brief get the recurrence coefficient \f$b_n\f$
          *
          * @param n get Nth coefficient
          */
         virtual DynamoFloat recurrenceB(const int n) const = 0;

         /**
          * @brief get the recurrence coefficient \f$c_n\f$
          *
          * @param n get Nth coefficient
          */
         virtual DynamoFloat recurrenceC(const int n) const = 0;

         /**
          * @brief Set the grid points
          */
         Array& rGrid();

         /**
          * @brief Set the polynomial values
          */
         Matrix& rPoly();

         /**
          * @brief Set the i-th derivative values
          *
          * @param i    order of the derivative
          */
         Matrix& rDiff(const int i);

         /**
          * @brief Initialise the polynomials from scratch
          */
         virtual void initFull() = 0;

         /**
          * @brief Initialise the polynomials with provided grid
          */
         virtual void initPartial() = 0;

         /**
          * @brief Compute the polynomials
          */
         virtual void computePolynomials() = 0;

         /**
          * @brief Compute the Derivatives
          */
         virtual void computeDerivatives() = 0;

      private:

         /**
          * @brief Number of grid points.
          */
         int mGridN;

         /**
          * @brief Number of polynomials in polynomial basis.
          */
         int mPolyN;

         /**
          * @brief Highest derivatives order.
          */
         int mMaxDiff;

         /**
          * @brief SmartArray of the grid points.
          */
         SmartArray    mpGrid;

         /**
          * @brief Matrix of the polynomial values on all grid points.
          */
         Matrix   mPoly;

         /**
          * @brief Vector of Matrices of the derivatives values on all grid points.
          */
         std::vector<Matrix>   mDiff;

         /**
          * @brief Create vector of Matrices for the derivatives
          */
         void initDiff();

   };

   // Public inline functions
   //
   inline int Polynomial::gridN() const
   {
      return this->mGridN;
   }

   inline int Polynomial::polyN() const
   {
      return this->mPolyN;
   }

   inline const Array& Polynomial::grid() const
   {
      return *this->mpGrid;
   }

   inline const SmartArray& Polynomial::pGrid() const
   {
      return this->mpGrid;
   }

   inline const Matrix& Polynomial::poly() const
   {
      return this->mPoly;
   }

   inline const Matrix& Polynomial::diff(const int i) const
   {
      return this->mDiff.at(i-1);
   }

   // Protected inline functions
   //

   inline int Polynomial::maxDiff()
   {
      return this->mMaxDiff;
   }

   inline Array& Polynomial::rGrid()
   {
      return *this->mpGrid;
   }

   inline Matrix& Polynomial::rPoly()
   {
      return this->mPoly;
   }

   inline Matrix& Polynomial::rDiff(const int i)
   {
      return this->mDiff.at(i-1);
   }

}

#endif // POLYNOMIAL_HPP
