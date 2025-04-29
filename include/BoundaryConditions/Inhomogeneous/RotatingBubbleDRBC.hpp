/** \file RotatingBubbleDRBC.hpp
 *  \brief Implementation of the second rotating bubble boundary condition
 */

#ifndef ROTATINGBUBBLEDRBC_HPP
#define ROTATINGBUBBLEDRBC_HPP

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
    *  \brief Implementation of the second rotating bubble boundary condition
    */
   class RotatingBubbleDRBC: public BoundaryCondition
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef SimulationConfig::NumericalScheme::RadialBasisType    BasisType;

         /**
          * @brief Constructor
          *
          * @param basis   Radial basis
          */
         RotatingBubbleDRBC(const BasisType &basis);

         /**
          * @brief Destructor
          */
         virtual ~RotatingBubbleDRBC() {};

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

#endif // ROTATINGBUBBLEDRBC_HPP
