/** \file GridGeneratorBase.hpp
 *  \brief Base of the grid generators implementations
 */

#ifndef GRIDGENERATORBASE_HPP
#define GRIDGENERATORBASE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the grid generators implementations
    */
   class GridGeneratorBase
   {
      public:
         /**
          * @brief Generate a Radial grid (\f$r\in [0, 1] \f$)
          */
         static const int RADIAL_GRID = 1;

         /**
          * @brief Generate a basic x grid (\f$ x\in [-1, 1] \f$)
          */
         static const int X_GRID = 2;

         /**
          * @brief Constructor
          *
          * @param gridN Number of grid points
          * @param type Type of the grid
          * @param isCSCS Generate a CSCS visualisation grid?
          */
         GridGeneratorBase(const int gridN, const int type, bool isCSCS = false);

         /**
          * @brief Destructor
          */
         virtual ~GridGeneratorBase() {};

         /**
          * @brief Get grid values
          */
         const Array& grid() const;

         /**
          * @brief Get quadrature weights 
          */
         const Array& weights() const;
         
      protected:
         /**
          * @brief Get size of the grid
          */
         int gridN() const;

         /**
          * @brief Convert \f$x\in [-1, 1]\f$ to radial \f$r \in [0, 1]\f$ grid
          */
         void convertToRadial();

         /**
          * @brief Sort grid (and weights) in increasing ordering
          */
         void sortGrid();

         /**
          * @brief Storage for the grid values
          */
         Array mGrid;

         /**
          * @brief Storage for the weights
          */
         Array mWeights;

         /**
          * @brief Grid type
          */
         int mType;

         /**
          * @brief Is a CSCS output type?
          */
         bool mIsCSCS;

      private:
   };

   inline int GridGeneratorBase::gridN() const
   {
      return this->mGrid.size();
   }

   inline const Array& GridGeneratorBase::grid() const
   {
      return this->mGrid;
   }
 
   inline const Array& GridGeneratorBase::weights() const
   {
      return this->mWeights;
   }

}

#endif // GRIDGENERATORBASE_HPP
