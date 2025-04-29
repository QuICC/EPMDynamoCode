/** \file MathConstants.hpp
 *  \brief Definition of some useful math constants
 */

#ifndef MATHCONSTANTS_HPP
#define MATHCONSTANTS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * @brief Contains some static useful math constants
    */
   class MathConstants
   {
      public:
         /**
         * @brief Empty Destructor
         */
         virtual ~MathConstants() {};

         /**
          * @brief The constant \f$\pi\f$
          */
         static const EPMFloat PI;

         /**
          * @brief Pure imaginary value I
          */
         static const EPMComplex cI;
         
      protected:

      private:
         /**
         * @brief Empty constructor
         */
         MathConstants() {};
   };

}

#endif // MATHCONSTANTS_HPP
