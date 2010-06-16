/** \file FDSHDegreeScalar.hpp
 *  \brief Implementation of the degree ordered FDSH data type
 */

#ifndef FDSHDEGREESCALAR_HPP
#define FDSHDEGREESCALAR_HPP

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"
#include "GeneralScalars/FDSHDegreeScalarBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the degree ordered FDSH data type
    */
   class FDSHDegreeScalar: public FDSHDegreeScalarBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param pTrunc Truncation information
          */
         FDSHDegreeScalar(const SmartTruncation pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~FDSHDegreeScalar() {};

         /**
          * @brief Get the coefficients stored by harmonic degree l
          *
          * This accesses the values with a fast access to the radial values 
          *
          * @param l   Requested  LShell degree l index
          */
         const MatrixZ&    lshell(const int l) const;

         /**
          * @brief Set the coefficients stored by harmonic degree l
          *
          * This accesses the values with a fast access to the radial values 
          *
          * @param l   Requested  LShell degree l index
          */
         MatrixZ&    rLShell(const int l);
         
      protected:

      private:
         /**
          * @brief Storage for the values in LShell ordering
          */
         std::vector<MatrixZ>  mLShells;

         /**
          * @brief Initialise the lshells
          */
         void initLShells();
   };

   inline const MatrixZ&   FDSHDegreeScalar::lshell(const int l) const
   {
      return this->mLShells.at(l);
   }

   inline MatrixZ&   FDSHDegreeScalar::rLShell(const int l)
   {
      return this->mLShells.at(l);
   }

}

#endif // FDSHDEGREESCALAR_HPP
