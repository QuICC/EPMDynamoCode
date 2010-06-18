/** \file RTPBase.hpp
 *  \brief Base for all the RTP space field implementation (scalar and vector)
 */

#ifndef RTPBASE_HPP
#define RTPBASE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"

namespace EPMDynamo {

   // Forward declaration of RTPField
   class RTPField;

   /**
    * @brief Basic structure required by any RTP field
    *
    * This provides some common methods required for all RTP fields
    */
   class RTPBase
   {
      public:
         /**
         * @brief Stores the pointer to a truncation object
         *
         * @param pTrunc truncation information
         */
         RTPBase(const SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~RTPBase() {};

         /**
          * @brief Get the number of spherical shells
          */
         int   nR() const;

         /**
          * @brief Get the number of theta coordinates
          *
          * @param n index of the concerned radial grid
          */
         int   nTh(const int n) const;

         /**
          * @brief Get the number of phi coordinates
          */
         int   nPh() const;

         /**
          * @brief Get radial grid point
          *
          * @param n Radial index
          */
         EPMFloat radGrid(const int n) const;

         /**
          * @brief Array of radial grid points
          */
         const Array&   radGrid() const;

         /**
          * @brief Get array of cos(theta) values
          *
          * @param j Theta index
          * @param n Radial index
          */
         EPMFloat cTh(const int j, const int n) const;

         /**
          * @brief Get array of sin(theta) values
          *
          * @param j Theta index
          * @param n Radial index
          */
         EPMFloat sTh(const int j, const int n) const;

         /**
          * @brief Get the Spectral Truncation pointer
          */
         SmartTruncation trunc() const;
         
      protected:
         /**
          * @brief Spectral truncation information
          */
         SmartTruncation   mpTrunc;

      private:
   };

   inline int   RTPBase::nR() const
   {
      return this->mpTrunc->local()->rtp()->nR();
   }

   inline int   RTPBase::nTh(const int n) const
   {
      return this->mpTrunc->local()->rtp()->nTh(n);
   }

   inline int   RTPBase::nPh() const
   {
      return this->mpTrunc->sim()->hoz()->nPh();
   }

   inline EPMFloat RTPBase::radGrid(const int n) const
   {
      return this->mpTrunc->local()->rtp()->radGrid(n);
   }

   inline const Array& RTPBase::radGrid() const
   {
      return this->mpTrunc->local()->rtp()->radGrid();
   }

   inline EPMFloat RTPBase::cTh(const int j, const int n) const
   {
      return this->mpTrunc->local()->rtp()->cTh(j, n);
   }

   inline EPMFloat RTPBase::sTh(const int j, const int n) const
   {
      return this->mpTrunc->local()->rtp()->sTh(j, n);
   }

   inline SmartTruncation RTPBase::trunc() const
   {
      return this->mpTrunc;
   }

}

#endif // RTPBASE_HPP
