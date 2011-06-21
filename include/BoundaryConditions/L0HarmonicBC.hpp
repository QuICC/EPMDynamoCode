/** \file L0HarmonicBC.hpp
 *  \brief General implementation of l=0 harmonic boundray condition
 */

#ifndef L0HARMONICBC_HPP
#define L0HARMONICBC_HPP

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "BoundaryConditions/BoundaryCondition.hpp"

namespace EPMDynamo {

   /**
    * @brief General description of a L0HarmonicBC
    *
    * This class defines the L0HarmonicBC concept independently of the used basis
    */
   class L0HarmonicBC: public BoundaryCondition
   {
      public:
         /**
          * @brief Constructor
          *
          * @param nL Number of harmonic degrees
          * @param nX Number of BC values
          * @param hasL0 CPU computes l=0 mode
          */
         L0HarmonicBC(const int nL, const int nX, const bool hasL0);

         /**
          * @brief Destructor
          */
         virtual ~L0HarmonicBC() {};

         /**
          * @brief Get boundary condition value RHS, real part
          *
          * @param l Harmonic degreee l
          * @param m Harmonic order m 
          */
         EPMFloat getRealRHSBC(const int l, const int m) const;

         /**
          * @brief Get boundary condition value RHS, imaginary part
          *
          * @param l Harmonic degreee l
          * @param m Harmonic order m 
          */
         EPMFloat getImagRHSBC(const int l, const int m) const;

         /**
          * @brief Get boundary condition value RHS, complex value
          *
          * @param l Harmonic degreee l
          * @param m Harmonic order m 
          */
         EPMComplex getRHSBC(const int l, const int m) const;
         
      protected:
         /**
          * @brief Flag used in parallel version to check if l=0 is stored locally
          */
         const bool mHasL0;

         /**
          * @brief Value of the l=0 mode
          */
         EPMFloat mValue;

      private:
   };

   inline EPMFloat L0HarmonicBC::getRealRHSBC(const int l, const int m) const
   {
      if(!this->mHasL0 || l != 0)
      {
         return 0.0;
      } else
      {
         return this->mValue;
      }
   }

   inline EPMFloat L0HarmonicBC::getImagRHSBC(const int l, const int m) const
   {
      return 0.0;
   }

   inline EPMComplex L0HarmonicBC::getRHSBC(const int l, const int m) const
   {
      if(!this->mHasL0 || l != 0)
      {
         return EPMComplex(0.0);
      } else
      {
         return EPMComplex(this->mValue);
      }
   }
}

#endif // L0HARMONICBC_HPP
