/** \file HorizontalDomain.hpp
 *  \brief Implementation of the simulation wide horizontal domain information (truncation + grids)
 */

#ifndef HORIZONTALDOMAIN_HPP
#define HORIZONTALDOMAIN_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Simulation/HorizontalTruncation.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the simulation wide horizontal domain information (truncation + grids)
    */
   class HorizontalDomain: public HorizontalTruncation
   {
      public:
         /**
          * @brief Constructor
          *
          * @param maxL  Highest harmonic degree
          * @param maxM  Highest harmonic order
          * @param mp    Longitudinal symmetry  assumption
          */
         HorizontalDomain(const int maxL, const int maxM, const int mp);

         /**
          * @brief Destructor
          */
         virtual ~HorizontalDomain() {};

         /**
          * @brief Get the \f$\theta\f$ grid
          */
         const Array&  thGrid() const;

         /**
          * @brief Get the \f$\cos\theta\f$ grid
          */
         const Array&  cosTheta() const;

         /**
          * @brief Get the \f$\sin\theta\f$ grid
          */
         const Array&  sinTheta() const;

         /**
          * @brief Get the \f$\varphi\f$ grid
          */
         const Array&  phGrid() const;

         /**
          * @brief Set the theta grid, cos theta and sin theta
          *
          * \param theta \f$\theta\f$ grid
          * \param cTh \f$\cos(\theta)\f$ grid
          * \param sTh \f$\sin(\theta)\f$ grid
          */
         void setThetaGrid(SmartArray theta, SmartArray cTh, SmartArray sTh);
         
      protected:

      private:
         /**
          * @brief Smart pointer to the \f$\theta\f$ grid
          */
         SmartArray  mpThGrid;

         /**
          * @brief Smart pointer to the \f$\cos\theta\f$ grid
          */
         SmartArray  mpCosTh;

         /**
          * @brief Smart pointer to the \f$\sin\theta\f$ grid
          */
         SmartArray  mpSinTh;

         /**
          * @brief Smart pointer to the \f$\varphi\f$ grid
          */
         SmartArray  mpPhGrid;

         /**
          * @brief Initialise the \f$\varphi\f$ grid from given truncation information
          */
         void initPhiGrid();
   };

   inline const Array& HorizontalDomain::thGrid() const
   {
      return *this->mpThGrid;
   }

   inline const Array& HorizontalDomain::cosTheta() const
   {
      return *this->mpCosTh;
   }

   inline const Array& HorizontalDomain::sinTheta() const
   {
      return *this->mpSinTh;
   }

   inline const Array& HorizontalDomain::phGrid() const
   {
      return *this->mpPhGrid;
   }

   inline void HorizontalDomain::setThetaGrid(SmartArray theta, SmartArray cTh, SmartArray sTh)
   {
      this->mpThGrid = theta;
      this->mpCosTh = cTh;
      this->mpSinTh = sTh;
   }

   /// Typedef for a smart shared pointer of HorizontalDomain
   typedef EPMSHARED_PTR<HorizontalDomain> SmartSimHozDomain;

}

#endif // HORIZONTALDOMAIN_HPP
