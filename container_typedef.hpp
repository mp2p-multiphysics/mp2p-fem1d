#ifndef CONTAINER_TYPEDEF
#define CONTAINER_TYPEDEF
#include <unordered_map>
#include <vector>

// typedef of vectors of numbers
typedef std::vector<double> VectorDouble;
typedef std::vector<int> VectorInt;

// typedef of unordered maps
typedef std::unordered_map<int, int> MapIntInt;

// nested vectors for domain integrals
typedef std::vector<double> Vector1D;
typedef std::vector<Vector1D> Vector2D;
typedef std::vector<Vector2D> Vector3D;
typedef std::vector<Vector3D> Vector4D;

// nested vectors and maps for boundary integrals
typedef std::unordered_map<int, std::unordered_map<int, double>> MapVector2D;
typedef std::unordered_map<int, std::unordered_map<int, Vector1D>> MapVector3D;
typedef std::unordered_map<int, std::unordered_map<int, Vector2D>> MapVector4D;

#endif
