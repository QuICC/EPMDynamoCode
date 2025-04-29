/** \file JacobiGrid.hpp
 *  \brief Implementation of a Jacobi grid generator
 */

#ifndef JACOBIGRID_HPP
#define JACOBIGRID_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "GridGenerators/GridGeneratorBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a Jacobi grid generator
    */
   class JacobiGrid: public GridGeneratorBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param gridN Size of the grid
          * @param type Type of the grid
          * @param isCSCS Generate a CSCS visualisation grid?
          */
         JacobiGrid(const int gridN, const int type = GridGeneratorBase::X_GRID, bool isCSCS = false);

         /**
          * @brief Destructor
          */
         virtual ~JacobiGrid() {};

         /**
          * @brief Recompute grid with a new set of alpha, beta parameters
          */
         void recompute(EPMFloat alpha, EPMFloat beta);
         
      protected:

      private:
         /**
          * @brief Alpha parameter
          */
         EPMFloat mAlpha;

         /**
          * @brief Beta parameter
          */
         EPMFloat mBeta;

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
          */
         void addCSCSRadial(int &left, int &right);
   };

}

#endif // JACOBIGRID_HPP
