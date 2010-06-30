/** \file InsulatingBC.hpp
 *  \brief Implementation of the magnetic insulating boundary condition
 */

#ifndef INSULATINGBC_HPP
#define INSULATINGBC_HPP

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
    * @brief This class defines the insulating homogeneous boundary condition
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class InsulatingBC: public HomogeneousBC
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::RadialBasisType    BasisType;

         /**
          * @brief Constructor
          *
          * @param basis   Radial basis
          */
         InsulatingBC(const BasisType &basis);

         /**
          * @brief Destructor
          */
         virtual ~InsulatingBC() {};

      protected:

         /**
          * @brief Fill in values for the boundaries
          *
          * @param basis   Radial basis
          */
         void fillLHSBCValues(const BasisType &basis);

      private:
   };

   template <typename TSimType> InsulatingBC<TSimType>::InsulatingBC(const typename InsulatingBC<TSimType>::BasisType &basis)
      : HomogeneousBC(basis.basisN(), basis.polyN())
   {
      // Fill Operator BC values
      this->fillLHSBCValues(basis);
   }

   template <typename TSimType> void InsulatingBC<TSimType>::fillLHSBCValues(const typename InsulatingBC<TSimType>::BasisType &basis)
   {
      int nL = this->nL();
      for(int l = 0; l < nL; ++l)
      {
         this->rLHSBCValues(l) = basis.at(l).bdiff(1) + static_cast<EPMFloat>(basis.at(l).l()+1)*basis.at(l).bpoly();
      }
   }

}

#endif // INSULATINGBC_HPP
