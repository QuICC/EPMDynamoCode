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
          * @brief Get the radial depend spherical geometry factor $\sqrt(l(l+1))$
          */
         const Array&  radSll() const;

         /**
          * @brief Set the radial grid
          *
          * \param rad Radial grid
          */
         void setRadGrid(SmartArray rad);

         /**
          * @brief Set the radial spherical geometric factor to stringentest value
          *
          * \param maxL Maximum harmonic degree
          */
         void initRadSll(int maxL);

         /**
          * @brief Set specific values for the radial spherical geometric factor
          */
         void setRadSll();
         
      protected:

      private:
         /**
          * @brief Smart pointer to the radial grid
          */
         SmartArray  mpRGrid;

         /**
          * @brief Smart pointer to the spherical geometric factor
          */
         SmartArray  mpRSll;
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
      this->mpRGrid = rad;
   }

   inline void RadialDomain::initRadSll(int maxL)
   {
      // Fill array with maximal value
      this->mpRSll->setConstant(std::sqrt(static_cast<EPMFloat>(maxL*(maxL+1))));
   }

   inline void RadialDomain::setRadSll()
   {
      // Currently set it by hand for the currently used resolution ...
      (*this->mpRSll)(0) = std::sqrt(2.0);
      (*this->mpRSll)(1) = std::sqrt(2.0);
      (*this->mpRSll)(2) = std::sqrt(12.0);
      (*this->mpRSll)(3) = std::sqrt(20.0);
      (*this->mpRSll)(4) = std::sqrt(42.0);
      (*this->mpRSll)(5) = std::sqrt(72.0);
      (*this->mpRSll)(6) = std::sqrt(110.0);
      (*this->mpRSll)(7) = std::sqrt(182.0);
      (*this->mpRSll)(8) = std::sqrt(272.0);
      (*this->mpRSll)(10) = std::sqrt(420.0);
      (*this->mpRSll)(11) = std::sqrt(552.0);
      (*this->mpRSll)(12) = std::sqrt(702.0);
   }

   /// Typedef for a smart simulation wide radial domain information 
   typedef EPMSHARED_PTR<RadialDomain> SmartSimRadDomain;

}

#endif // RADIALDOMAIN_HPP
