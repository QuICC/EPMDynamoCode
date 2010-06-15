/** \file DifferentialHeatingBC.hpp
 *  \brief Implemenation of the l=0 differential heating boundary condition
 */

#ifndef DIFFERENTIALHEATINGBC_HPP
#define DIFFERENTIALHEATINGBC_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "BoundaryConditions/L0HarmonicBC.hpp"

namespace EPMDynamo {

   /**
    * @brief This class defines the  L=0 harmonic differential heating boundary condition
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class DifferentialHeatingBC: public L0HarmonicBC
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::RadialBasisType    BasisType;

         /**
          * @brief Constructor
          *
          * @param bcValue Value at the boundary
          * @param basis   Radial basis
          */
         DifferentialHeatingBC(const DynamoFloat bcValue, const BasisType &basis);

         /**
          * @brief Destructor
          */
         virtual ~DifferentialHeatingBC() {};

      protected:

         /**
          * @brief Fill in values for the boundaries
          *
          * @param basis   Radial basis
          */
         void fillLHSBCValues(const BasisType &basis);

         /**
          * @brief Set Boundary value
          *
          * @param bcValue   Value at the boundary
          */
         void setBCValue(const DynamoFloat bcValue);

      private:
   };

   template <typename TSimType> DifferentialHeatingBC<TSimType>::DifferentialHeatingBC(const DynamoFloat bcValue, const typename DifferentialHeatingBC<TSimType>::BasisType &basis)
      : L0HarmonicBC(basis.basisN(), basis.polyN())
   {
      // Fill Operator BC values
      this->fillLHSBCValues(basis);

      // Set boundary value
      this->setBCValue(bcValue);

      // Check if boundary condition is homogeneous
      if(bcValue == 0.0)
      {
         this->mIsHomogeneous = true;
      }
   }

   template <typename TSimType> void DifferentialHeatingBC<TSimType>::fillLHSBCValues(const typename DifferentialHeatingBC<TSimType>::BasisType &basis)
   {
      int nL = this->nL();
      for(int l = 0; l < nL; ++l)
      {
         this->rLHSBCValues(l) = basis.at(l).bpoly();
      }
   }

   template <typename TSimType> void DifferentialHeatingBC<TSimType>::setBCValue(const DynamoFloat bcValue)
   {
      this->mValue = bcValue;
   }

}

#endif // DIFFERENTIALHEATINGBC_HPP
