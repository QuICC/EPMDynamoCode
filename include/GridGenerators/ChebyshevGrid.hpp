/** \file ChebyshevGrid.hpp
 *  \brief Implementation of a Chebyshev grid generator
 */

#ifndef CHEBYSHEVGRID_HPP
#define CHEBYSHEVGRID_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "GridGenerators/GridGeneratorBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a Chebyshev grid generator
    */
   class ChebyshevGrid: public GridGeneratorBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param gridN Size of the grid
          * @param type Type of the generated grid
          * @param isCSCS Generate a CSCS visualisation grid?
          */
         ChebyshevGrid(const int gridN, const int type = GridGeneratorBase::X_GRID, bool isCSCS = false);

         /**
          * @brief Destructor
          */
         virtual ~ChebyshevGrid() {};
         
      protected:

      private:
         /**
          * @brief Generate grid and weights
          */
         void initialise();

         /**
          * @brief Generate grid and weights
          *
          * @param left Left range offset for partial grid
          * @param right Right range offset for partial grid
          */
         void computeWeightedGrid(int left = 0, int right = 0);

         /**
          * @brief Generate the grid points
          *
          * @param left Left range offset for partial grid
          * @param right Right range offset for partial grid
          */
         void computeGrid(const int left, const int right);

         /**
          * @brief Generat the quadrature weights
          *
          * @param left Left range offset for partial grid
          * @param right Right range offset for partial grid
          */
         void computeWeights(const int left, const int right);

         /**
          * @brief Include additional points for CSCS radial grid
          *
          * @param left Left offset
          * @param right Right offset
          */
         void addCSCSRadial(int &left, int &right);
   };

}

#endif // CHEBYSHEVGRID_HPP
