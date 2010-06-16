/** \file FDSHShell.hpp
 *  \brief Implementation of a generic "FDSH spherical shell"
 */

#ifndef FDSHSHELL_HPP
#define FDSHSHELL_HPP

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"

namespace EPMDynamo {

   /**
    * @brief Abstraction for a "harmonic shell" in a FD Spherical harmonics decomposition
    *
    * This type of object is required to store the coefficients for the radius
    * independent operations on the harmomics modes.
    *
    * The best storage orders seems to be store the coefficients by array of degrees
    * per order. These values can't be stored in a matrix (or at least not without
    * wasting lots of memory) as they have different lengths.
    */
   class FDSHShell
   {
      public:
         /**
         * @brief Constructs the shell based on the given truncation information
         *
         * @param pTrunc Truncation information pointer
         */
         FDSHShell(SmartTruncation pTrunc);

         /**
         * @brief Simple empty destructor
         */
         virtual ~FDSHShell() {};

         /**
          * @brief Get the truncation information
          */
         const SmartTruncation trunc() const;

         /**
          * @brief   Get array of coefficient for given order
          *
          * @param   m  Order of the requested SH coefficients
          */
         const SHDegrees&   order(const int   m) const;

         /**
          * @brief   Return array of coefficient for given order
          *
          * @param   m  Order of the requested SH coefficients
          */
         SHDegrees&   rOrder(const int   m);
         
      protected:

      private:

         /**
          * @brief Pointer to truncation information
          */
         SmartTruncation  mpTrunc;

         /** 
          * @brief Storage for the modes. it's an array of degrees for each order.
          */
         std::vector<SHDegrees>    mOrders;

         /**
          * @brief Initialise the orders
          */
         void initOrders();
   };

   inline const SmartTruncation FDSHShell::trunc() const
   {
      return this->mpTrunc;
   }

   inline const SHDegrees& FDSHShell::order(const int m) const
   {
      return this->mOrders.at(m);
   }

   inline SHDegrees& FDSHShell::rOrder(const int m)
   {
      return mOrders.at(m);
   }

}

#endif // FDSHSHELL_HPP
