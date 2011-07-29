/** \file TimeDependentBC.hpp
 *  \brief General implementation for a time dependent boundary condition
 */

#ifndef TIMEDEPENDENTBC_HPP
#define TIMEDEPENDENTBC_HPP

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "BoundaryConditions/BoundaryCondition.hpp"
#include "Timestepping/TimestepParameters.hpp"

namespace EPMDynamo {

   /**
    * \brief General implementation for a time dependent boundary condition
    */
   class TimeDependentBC: public BoundaryCondition
   {
      public:
         /**
          * @brief Constructor
          *
          * @param nL Number of harmonic degrees
          * @param nX Number of BC values
          */
         TimeDependentBC(const int nL, const int nX, const TimestepParameters &tsParams);

         /**
          * @brief Destructor
          */
         virtual ~TimeDependentBC() {};

         /**
          * @brief Get boundary condition value RHS, real part
          *
          * @param l Harmonic degreee l
          * @param m Harmonic order m 
          */
         EPMFloat getRealRHSBC(const int l, const int m) const = 0;

         /**
          * @brief Get boundary condition value RHS, imaginary part
          *
          * @param l Harmonic degreee l
          * @param m Harmonic order m 
          */
         EPMFloat getImagRHSBC(const int l, const int m) const = 0;

         /**
          * @brief Get boundary condition value RHS, complex value
          *
          * @param l Harmonic degreee l
          * @param m Harmonic order m 
          */
         EPMComplex getRHSBC(const int l, const int m) const = 0;
         
      protected:
         /**
          * @brief Function providing the time dependent real part factor
          */
         virtual EPMFloat evolvingRealFactor() const = 0;

         /**
          * @brief Function providing the time dependent imaginary part factor
          */
         virtual EPMFloat evolvingImagFactor() const = 0;

         /**
          * @brief Reference to a TimestepParameters object to get access to time and timestep
          */
         const TimestepParameters&  mrTSParams;

      private:
   };
}

#endif // TIMEDEPENDENTBC_HPP
