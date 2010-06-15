/** \file StressFreeTorBC.hpp
 *  \brief Implemenation of the stress free velocity boundary condition on the toroidal scalar
 */

#ifndef STRESSFREETORBC_HPP
#define STRESSFREETORBC_HPP

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
    * @brief This class defines the stress free homogeneous boundary condition on a
    *    toroidal scalar
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class StressFreeTorBC: public HomogeneousBC
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::RadialBasisType    BasisType;

         /**
          * @brief Constructor
          *
          * @param basis   Radial basis
          */
         StressFreeTorBC(const BasisType &basis);

         /**
          * @brief Destructor
          */
         virtual ~StressFreeTorBC() {};

      protected:

         /**
          * @brief Fill in values for the boundaries
          *
          * @param basis   Radial basis
          */
         void fillLHSBCValues(const BasisType &basis);

      private:
   };

   template <typename TSimType> StressFreeTorBC<TSimType>::StressFreeTorBC(const typename StressFreeTorBC<TSimType>::BasisType &basis)
      : HomogeneousBC(basis.basisN(), basis.polyN())
   {
      // Fill Operator BC values
      fillLHSBCValues(basis);
   }

   template <typename TSimType> void StressFreeTorBC<TSimType>::fillLHSBCValues(const typename StressFreeTorBC<TSimType>::BasisType &basis)
   {
      int nL = this->nL();
      for(int l = 0; l < nL; ++l)
      {
         rLHSBCValues(l) = basis.at(l).bdiff(1) - basis.at(l).bpoly();
      }
   }

}

#endif // STRESSFREETORBC_HPP
