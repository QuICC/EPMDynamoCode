/** \file WorlandChebyshevTraits.hpp
 *  \brief Traits to use the Chebyshev grid Worland polynomials
 */

#ifndef WORLANDCHEBYSHEVTRAITS_HPP
#define WORLANDCHEBYSHEVTRAITS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "GridGenerators/ChebyshevGrid.hpp"

namespace EPMDynamo {

   /**
    * @brief Traits to use the Chebyshev grid Worland polynomials
    */
   class WorlandChebyshevTraits
   {
      public:
         /// Typedef for setting the radial grid
         typedef ChebyshevGrid RadialGrid;

         /**
          * @brief Alpha parameter of the Jacobi polynomials for a Worland Chebyshev quadrature rule
          */
         static const DynamoFloat ALPHA;

         /**
          * @brief beta parameter of the Jacobi polynomials for a Worland Chebyshev quadrature rule
          *
          * This is only part of the parameter, the actual parameter is \f$\beta = l + BETA\f$
          */
         static const DynamoFloat BETA;
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         WorlandChebyshevTraits();

         /**
          * @brief Destructor
          */
         virtual ~WorlandChebyshevTraits() {};
   };

}

#endif // WORLANDCHEBYSHEVTRAITS_HPP
