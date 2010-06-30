/** \file DRadialBC.hpp
 *  \brief Implementation of the homogeneous first radial derivative boundary condition
 */

#ifndef DRADIALBC_HPP
#define DRADIALBC_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "BoundaryConditions/HomogeneousBC.hpp"

namespace EPMDynamo {

   /**
    * @brief This class defines the homogeneous boundary condition on the first radial
    *    derivative
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class DRadialBC: public HomogeneousBC
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::RadialBasisType    BasisType;

         /**
          * @brief Constructor
          *
          * @param basis   Radial basis
          */
         DRadialBC(const BasisType &basis);

         /**
          * @brief Destructor
          */
         virtual ~DRadialBC() {};

      protected:

         /**
          * @brief Fill in values for the boundaries
          *
          * @param basis   Radial basis
          */
         void fillLHSBCValues(const BasisType &basis);

      private:
   };

   template <typename TSimType> DRadialBC<TSimType>::DRadialBC(const typename DRadialBC<TSimType>::BasisType &basis)
      : HomogeneousBC(basis.basisN(), basis.polyN())
   {
      // Fill Operator BC values
      this->fillLHSBCValues(basis);
   }

   template <typename TSimType> void DRadialBC<TSimType>::fillLHSBCValues(const typename DRadialBC<TSimType>::BasisType &basis)
   {
      int nL = this->nL();
      for(int l = 0; l < nL; ++l)
      {
         this->rLHSBCValues(l) = basis.at(l).bdiff(1);
      }
   }
}

#endif // DRADIALBC_HPP
