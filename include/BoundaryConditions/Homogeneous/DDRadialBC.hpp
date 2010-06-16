/** \file DDRadialBC.hpp
 *  \brief Implementation of the homogeneous second radial derivative boundary condition
 */

#ifndef DDRADIALBC_HPP
#define DDRADIALBC_HPP

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
    * @brief This class defines the homogeneous boundary condition on the second radial
    *    derivative
    *
    * \tparam TSimTypeType Type of the simulation
    */
   template <typename TSimType> class DDRadialBC: public HomogeneousBC
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::RadialBasisType    BasisType;

         /**
          * @brief Constructor
          *
          * @param basis   Radial basis
          */
         DDRadialBC(const BasisType &basis);

         /**
          * @brief Destructor
          */
         virtual ~DDRadialBC() {};

      protected:

         /**
          * @brief Fill in values for the boundaries
          *
          * @param basis   Radial basis
          */
         void fillLHSBCValues(const BasisType &basis);

      private:
   };

   template <typename TSimType> DDRadialBC<TSimType>::DDRadialBC(const typename DDRadialBC<TSimType>::BasisType &basis)
      : HomogeneousBC(basis.basisN(), basis.polyN())
   {
      // Fill Operator BC values
      fillLHSBCValues(basis);
   }

   template <typename TSimType> void DDRadialBC<TSimType>::fillLHSBCValues(const typename DDRadialBC<TSimType>::BasisType &basis)
   {
      int nL = this->nL();
      for(int l = 0; l < nL; ++l)
      {
         rLHSBCValues(l) = basis.at(l).bdiff(2);
      }
   }
}

#endif // DDRADIALBC_HPP
