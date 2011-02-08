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

         /**
          * @brief Get the radial depend spherical geometry factor $\sqrt(l(l+1))$
          */
         const Array&  radSll() const;

         /**
          * @brief Set the radial spherical geometric factor to stringentest value
          *
          * \param maxL Maximum harmonic degree
          */
         void initRadSll(int maxL);
         
      protected:

      private:
         /**
          * @brief Local storage for the highest L
          */
         int   mMaxL;

         /**
          * @brief Smart pointer to the radial grid
          */
         SmartArray  mpRGrid;

         /**
          * @brief Smart pointer to the spherical geometric factor
          */
         SmartArray  mpRSll;

         /**
          * @brief Set specific values for the radial spherical geometric factor
          */
         void setRadSll();

         /**
          * @brief Compute the alpha parameter of the Jacobi polynomials
          */
         EPMFloat alpha() const;

         /**
          * @brief Compute the beta parameter of the Jacobi polynomials
          */
         EPMFloat beta(int l) const;
   };

   inline const Array& RadialDomain::radGrid() const
   {
      return *this->mpRGrid;
   }

   inline const Array& RadialDomain::radSll() const
   {
      return *this->mpRSll;
   }

   inline void RadialDomain::setRadGrid(SmartArray rad)
   {
      // Set the radial grid
      this->mpRGrid = rad;

      // Set the sqrt(l(l+1)) factors related to the grid
      this->setRadSll();
   }

   inline void RadialDomain::initRadSll(int maxL)
   {
      // Store the maximul harmonic degree
      this->mMaxL = maxL;

      // Fill array with maximal value
      this->mpRSll->setConstant(std::sqrt(static_cast<EPMFloat>(maxL*(maxL+1))));
   }

   /// Typedef for a smart simulation wide radial domain information 
   typedef EPMSHARED_PTR<RadialDomain> SmartSimRadDomain;

}

#endif // RADIALDOMAIN_HPP
