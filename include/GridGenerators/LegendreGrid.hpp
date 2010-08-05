/** \file LegendreGrid.hpp
 *  \brief Implementation of a Legendre grid generator
 */

#ifndef LEGENDREGRID_HPP
#define LEGENDREGRID_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "GridGenerators/GridGeneratorBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a Legendre grid generator
    */
   class LegendreGrid: public GridGeneratorBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param gridN Size of the grid
          * @param type Type of the grid
          * @param isCSCS Generate a CSCS visualisation grid?
          */
         LegendreGrid(const int gridN, const int type = GridGeneratorBase::X_GRID, bool isCSCS = false);

         /**
          * @brief Destructor
          */
         virtual ~LegendreGrid() {};
         
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
          * @brief Compute the quadrature matrix
          *
          * @param diag Diagonal of the quadrature matrix
          * @param subdiag Subdiagonal of the quadrature matrix
          */
         void computeQuadratureMatrix(Array& diag, Array& subdiag);

         /**
          * @brief Do modifications for CSCS output
          *
          * @param left Left offset
          * @param right Right offset
          */
         void addCSCS(int &left, int &right);

         /**
          * @brief Include additional points for CSCS theta grid
          *
          * @param left Left offset
          * @param right Right offset
          */
         void addCSCSTheta(int &left, int &right);

         /**
          * @brief Include additional points for CSCS radial grid
          *
          * @param left Left offset
          * @param right Right offset
          *
          * \epmTodo The Legendre CSCS grid seems to only add r=1 and not r=0!
          */
         void addCSCSRadial(int &left, int &right);
   };

}

#endif // LEGENDREGRID_HPP
