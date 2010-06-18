/** \file FDSHOrderScalar.hpp
 *  \brief Implementation of the harmonic order ordered FDSH data type
 */

#ifndef FDSHORDERSCALAR_HPP
#define FDSHORDERSCALAR_HPP

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"
#include "GeneralScalars/FDSHOrderScalarBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the harmonic order ordered FDSH data type
    */
   class FDSHOrderScalar: public FDSHOrderScalarBase
   {
      public:
         /**
          * @brief Constructor
          */
         FDSHOrderScalar(const SmartTruncation pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~FDSHOrderScalar() {};

         /**
          * @brief Get the harmonic order "spherical FDSH shell"
          *
          * @param m Harmonic order m
          */
         const MatrixZ&    mshell(const int m) const;

         /**
          * @brief Set the harmonic order "spherical FDSH shell"
          *
          * @param m Harmonic order m
          */
         MatrixZ&    rMShell(const int m);

         /**
          * @brief Get coefficient for mode parameters l, m, and n
          *
          * This routine should only be used if nothing else is possible. It doesn't show anything about 
          * fast or slow looping directions. It is only here for a few special case where it allows a more
          * general implementation.
          */
         EPMComplex  coeff(const int l, const int m, const int n) const;

         /**
          * @brief Set coefficient for mode parameters l, m, and n
          *
          * This routine should only be used if nothing else is possible. It doesn't show anything about 
          * fast or slow looping directions. It is only here for a few special case where it allows a more
          * general implementation. This might be removed at a later stage.
          */
         EPMComplex&  rCoeff(const int l, const int m, const int n);
         
      protected:

      private:
         /**
          * @brief Storage for the values in MShell ordering
          */
         std::vector<MatrixZ>  mMShells;

         /**
          * @brief Initialise the mshells
          */
         void initMShells();
   };

   inline const MatrixZ&   FDSHOrderScalar::mshell(const int m) const
   {
      return this->mMShells.at(m);
   }

   inline MatrixZ&   FDSHOrderScalar::rMShell(const int m)
   {
      return this->mMShells.at(m);
   }

   inline EPMComplex  FDSHOrderScalar::coeff(const int l, const int m, const int n) const
   {
      return this->mMShells.at(m)(l, n);
   }

   inline EPMComplex&  FDSHOrderScalar::rCoeff(const int l, const int m, const int n)
   {
      return this->mMShells.at(m)(l, n);
   }

}

#endif // FDSHORDERSCALAR_HPP
