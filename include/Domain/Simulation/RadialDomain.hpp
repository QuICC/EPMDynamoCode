/** \file RadialDomain.hpp
 *  \brief Implementation of the simulation wide radial domain information (truncation + grids)
 */

#ifndef RADIALDOMAIN_HPP
#define RADIALDOMAIN_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Simulation/RadialTruncation.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the simulation wide radial domain information (truncation + grids)
    */
   class RadialDomain: public RadialTruncation
   {
      public:
         /**
          * @brief Constructor
          *
          * @param maxN Highest radial expansion order
          * @param nR Number of radial grid points
          */
         RadialDomain(const int maxN, const int nR);

         /**
          * @brief Destructor
          */
         virtual ~RadialDomain() {};

         /**
          * @brief Get the radial grid
          */
         const Array&  radGrid() const;

         /**
          * @brief Set the radial grid
          *
          * \param rad Radial grid
          */
         void setRadGrid(SmartArray rad);
         
      protected:

      private:
         /**
          * @brief Smart pointer to the radial grid
          */
         SmartArray  mpRGrid;
   };

   inline const Array& RadialDomain::radGrid() const
   {
      return *this->mpRGrid;
   }

   inline void RadialDomain::setRadGrid(SmartArray rad)
   {
      this->mpRGrid = rad;
   }

   /// Typedef for a smart simulation wide radial domain information 
   typedef EPMSHARED_PTR<RadialDomain> SmartSimRadDomain;

}

#endif // RADIALDOMAIN_HPP
