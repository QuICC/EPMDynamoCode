/** \file RTPDomain.hpp
 *  \brief Implementation of the remote (core) wide RTP space domain information
 */

#ifndef RTPDOMAIN_HPP
#define RTPDOMAIN_HPP

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
#include "Domain/Core/RTPTruncation.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the remote (core) wide RTP space domain information
    */
   class RTPDomain: public RTPTruncation
   {
      public:
         /**
          * @brief Constructor
          *
          * \param r0 Radial start index
          * \param nR Number of radial indexes
          * \param th0 Array of theta start indexes
          * \param nTh Array of number of theta values
          */
         RTPDomain(const int r0, const int nR, const ArrayI &th0, const ArrayI &nTh);

         /**
          * @brief Destructor
          */
         virtual ~RTPDomain() {};

         /**
          * @brief Get a radial grid point
          *
          * \param n Radial index
          */
         DynamoFloat radGrid(const int n) const;

         /**
          * @brief Get the radial grid
          */
         const Array& radGrid() const;

         /**
          * @brief Get the j \f$\cos\theta\f$ value (for radial index n)
          *
          * \param j Theta index
          * \param n Radial index 
          */
         DynamoFloat  cTh(const int j, const int n) const;

         /**
          * @brief Get the j \f$\sin\theta\f$ value (for radial index n)
          *
          * \param j Theta index
          * \param n Radial index 
          */
         DynamoFloat  sTh(const int j, const int n) const;

         /**
          * @brief Set the radial grid (isn't know at creation time)
          *
          * @param rad Radial grid
          */
         void setRadialGrid(SmartArray rad);

         /**
          * @brief Set cos theta and sin theta
          *
          * \param cTh \f$\cos(\theta)\f$ grid
          * \param sTh \f$\sin(\theta)\f$ grid
          */
         void setThetaGrid(SmartArray cTh, SmartArray sTh);
         
      protected:

      private:
         /**
          * @brief Smart pointer to the radial grid
          */
         SmartArray mpRadGrid;

         /**
          * @brief Smart pointer to the \f$\cos\theta\f$ grid
          */
         SmartArray  mpCosTh;

         /**
          * @brief Smart pointer to the \f$\sin\theta\f$ grid
          */
         SmartArray  mpSinTh;
   };

   inline DynamoFloat RTPDomain::radGrid(const int n) const
   {
      return (*this->mpRadGrid)(n);
   }

   inline const Array& RTPDomain::radGrid() const
   {
      return *this->mpRadGrid;
   }

   inline DynamoFloat RTPDomain::cTh(const int j, const int n) const
   {
      return (*this->mpCosTh)(this->th0(n) + j);
   }

   inline DynamoFloat RTPDomain::sTh(const int j, const int n) const
   {
      return (*this->mpSinTh)(this->th0(n) + j);
   }

   /// Typedef for a smart RTPDomain
   typedef EPMSHARED_PTR<RTPDomain> SmartCoreRTPDomain;

}

#endif // RTPDOMAIN_HPP
