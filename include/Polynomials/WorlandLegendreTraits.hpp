/** \file WorlandLegendreTraits.hpp
 *  \brief Traits to use the Legendre grid Worland polynomials
 */

#ifndef WORLANDLEGENDRETRAITS_HPP
#define WORLANDLEGENDRETRAITS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "GridGenerators/LegendreGrid.hpp"

namespace EPMDynamo {

   /**
    * @brief Traits to use the Legendre grid Worland polynomials
    */
   class WorlandLegendreTraits
   {
      public:
         /// Typedef for setting the radial grid
         typedef LegendreGrid RadialGrid;

         /**
          * @brief Alpha parameter of the Jacobi polynomials for a Worland Legendre quadrature rule
          */
         static const DynamoFloat ALPHA;

         /**
          * @brief beta parameter of the Jacobi polynomials for a Worland Legendre quadrature rule
          *
          * This is only part of the parameter, the actual parameter is \f$\beta = l + BETA\f$
          */
         static const DynamoFloat BETA;
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         WorlandLegendreTraits();

         /**
          * @brief Destructor
          */
         virtual ~WorlandLegendreTraits() {};
   };

}

#endif // WORLANDLEGENDRETRAITS_HPP
