/** \file RadialTruncation.hpp
 *  \brief Implementation of the simulation radial truncation information
 */

#ifndef RADIALTRUNCATION_HPP
#define RADIALTRUNCATION_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the simulation radial truncation information
    */
   class RadialTruncation
   {
      public:
         /**
          * @brief Constructor
          *
          * @param maxN Highest radial expansion order
          * @param nR Number of radial grid points
          */
         RadialTruncation(const int maxN, const int nR);

         /**
          * @brief Destructor
          */
         virtual ~RadialTruncation() {};

         /**
          * @brief Get the spectral maximum radial index
          */
         int maxN() const;

         /**
          * @brief Get the number of radial spectral coefficients
          */
         int nN() const;

         /**
          * @brief Get the humber of radial grid points
          */
         int nR() const;
         
      protected:

      private:
         /**
          * @brief Maximum spectral radial index
          */
         int mMaxN;

         /**
          * @brief Number of spectral coefficients
          */
         int mNn;

         /**
          * @brief Size of real space radial grid
          */
         int mNr;

         /**
          * @brief Correct grid sizes if grid type requires it (ex. CSCS)
          */
         void correctGridSizes();
   };

   inline int RadialTruncation::maxN() const
   {
      return this->mMaxN;
   }

   inline int RadialTruncation::nN() const
   {
      return this->mNn;
   }

   inline int RadialTruncation::nR() const
   {
      return this->mNr;
   }

   /// Typedef for a smart simulation wid radial truncation information
   typedef EPMSHARED_PTR<RadialTruncation> SmartSimRadTrunc;

}

#endif // RADIALTRUNCATION_HPP
