#ifndef CONTAINER_TYPEDEF
#define CONTAINER_TYPEDEF
#include <vector>

// typedef of vectors of numbers
typedef std::vector<double> VectorDouble;
typedef std::vector<int> VectorInt;

// nested vectors for integration
typedef std::vector<double> Vector1D;
typedef std::vector<Vector1D> Vector2D;
typedef std::vector<Vector2D> Vector3D;
typedef std::vector<Vector3D> Vector4D;

#endif
