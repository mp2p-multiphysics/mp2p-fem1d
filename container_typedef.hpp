#ifndef CONTAINER_TYPEDEF
#define CONTAINER_TYPEDEF
#include <unordered_map>
#include <vector>
#include "Eigen/Eigen"

namespace FEM1D
{

// nested int vectors
typedef std::vector<int> VectorInt;
typedef std::vector<VectorInt> VectorInt2D;

// nested double vectors
typedef std::vector<double> VectorDouble;
typedef std::vector<VectorDouble> VectorDouble2D;
typedef std::vector<VectorDouble2D> VectorDouble3D;
typedef std::vector<VectorDouble3D> VectorDouble4D;

// unordered maps
typedef std::unordered_map<int, int> MapIntInt;

// Eigen objects
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> EigenSparseMatrix;
typedef Eigen::Triplet<double> EigenTriplet;
typedef Eigen::VectorXd EigenVector;

}

#endif
