/** \file LatitudinalLibrationBC.hpp
 *  \brief Implementation of the time dependent toroidal latitudinal libration boundary condition
 */

#ifndef LATITUDINALLIBRATIONBC_HPP
#define LATITUDINALLIBRATIONBC_HPP

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
    * \brief Implementation of the time dependent toroidal latitudinal libration boundary condition
    */
   class LatitudinalLibrationBC: public TimeDependentBC
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef SimulationConfig::NumericalScheme::RadialBasisType    BasisType;

         /**
          * @brief Constructor
          *
          * @param basis   Radial basis
          */
         LatitudinalLibrationBC(const BasisType &basis, const TimestepParameters &tsParams);

         /**
          * @brief Destructor
          */
         virtual ~LatitudinalLibrationBC() {};

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
         EPMFloat mEpsilonx;

         /**
          * @brief Epsilon parameter
          */
         EPMFloat mEpsilony;
         
         /**
          * @brief Frequency parameter
          */
         EPMFloat mFrequency;
   };
}

#endif // LATITUDINALLIBRATIONBC_HPP
