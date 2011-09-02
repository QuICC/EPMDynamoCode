/** \file LongitudinalLibrationBC.hpp
 *  \brief Implementation of the time dependent toroidal longitudinal libration boundary condition
 */

#ifndef LONGITUDINALLIBRATIONBC_HPP
#define LONGITUDINALLIBRATIONBC_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"
#include "Config/NumericalSchemeInc.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "BoundaryConditions/TimeDependentBC.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the time dependent toroidal longitudinal libration boundary condition
    *
    * \epmBug Need to check the amplitude of the m=0 mode
    */
   class LongitudinalLibrationBC: public TimeDependentBC 
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef SimulationConfig::NumericalScheme::RadialBasisType    BasisType;

         /**
          * @brief Constructor
          *
          * @param basis   Radial basis
          */
         LongitudinalLibrationBC(const BasisType &basis, const TimestepParameters &tsParams);

         /**
          * @brief Destructor
          */
         virtual ~LongitudinalLibrationBC() {};

         /**
          * @brief Get boundary condition value RHS, real part
          *
          * @param l Harmonic degreee l
          * @param m Harmonic order m 
          */
         virtual EPMFloat getRealRHSBC(const int l, const int m) const;

         /**
          * @brief Get boundary condition value RHS, imaginary part
          *
          * @param l Harmonic degreee l
          * @param m Harmonic order m 
          */
         virtual EPMFloat getImagRHSBC(const int l, const int m) const;

         /**
          * @brief Get boundary condition value RHS, complex value
          *
          * @param l Harmonic degreee l
          * @param m Harmonic order m 
          */
         virtual EPMComplex getRHSBC(const int l, const int m) const;

      protected:
         /**
          * @brief Function providing the time dependent real part factor
          */
         virtual EPMFloat evolvingRealFactor() const;

         /**
          * @brief Function providing the time dependent imaginary part factor
          */
         virtual EPMFloat evolvingImagFactor() const;

         /**
          * @brief Setup the parametrisation
          */
         void setup();

         /**
          * @brief Fill in values for the boundaries
          *
          * @param basis   Radial basis
          */
         void fillLHSBCValues(const BasisType &basis);

      private:
         /**
          * @brief Epsilon parameter
          */
         EPMFloat mEpsilon;

         /**
          * @brief Frequency parameter
          */
         EPMFloat mFrequency;
   };
}

#endif // LONGITUDINALLIBRATIONBC_HPP
