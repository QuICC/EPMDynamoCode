/** \file PrecessionFrameBC.hpp
 *  \brief Implementation of toroidal BC for precession forcing
 */

#ifndef PRECESSIONFRAMEBC_HPP
#define PRECESSIONFRAMEBC_HPP

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
#include "BoundaryConditions/BoundaryCondition.hpp"

namespace EPMDynamo {

   /**
    *  \brief Implementation of toroidal BC for precession forcing
    *
    *  This boundary condition is not timedepent because the frame of reference is moved so as to be able to use
    *  time independent BCs.
    */
   class PrecessionFrameBC: public BoundaryCondition
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef SimulationConfig::NumericalScheme::RadialBasisType    BasisType;

         /**
          * @brief Constructor
          *
          * @param basis   Radial basis
          */
         PrecessionFrameBC(const BasisType &basis);

         /**
          * @brief Destructor
          */
         virtual ~PrecessionFrameBC() {};

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
          * @brief Fill in values for the boundaries
          *
          * @param basis   Radial basis
          */
         void fillLHSBCValues(const BasisType &basis);

      private:
   };
}

#endif // PRECESSIONFRAMEBC_HPP
