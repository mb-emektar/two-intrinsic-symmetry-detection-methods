#pragma once

#define _CRT_SECURE_NO_DEPRECATE
#include <inttypes.h>
#define HAVE_INT8_T
#include <cmath>
#include <xlocnum>
#define _USE_MATH_DEFINES // For M_PI

#include <Inventor/SbVec3f.h>
#include <iostream>
#include <vector>
#include <numeric> // Required for std::accumulate

namespace mesh {
	
#include <ostream>  // For std::ostream (used in operator<<)
#include <cassert>  // For assert (used in operator[])

	struct Vec3 {
		float x, y, z;

		// Constructors
		Vec3() : x(0), y(0), z(0) {};
		Vec3(float x1, float y1, float z1) : x(x1), y(y1), z(z1) {};

		// --- Member-wise Assignment Operators (modifies the vector) ---

		// Vector assignment addition
		Vec3& operator+=(const Vec3& rhs) {
			x += rhs.x;
			y += rhs.y;
			z += rhs.z;
			return *this;
		}

		// Vector assignment subtraction
		Vec3& operator-=(const Vec3& rhs) {
			x -= rhs.x;
			y -= rhs.y;
			z -= rhs.z;
			return *this;
		}

		// Scalar assignment multiplication
		Vec3& operator*=(float scalar) {
			x *= scalar;
			y *= scalar;
			z *= scalar;
			return *this;
		}

		// Scalar assignment division
		Vec3& operator/=(float scalar) {
			x /= scalar;
			y /= scalar;
			z /= scalar;
			return *this;
		}

		// --- Unary Operator ---

		// Unary negation
		Vec3 operator-() const {
			return { -x, -y, -z };
		}

		// --- Indexing ---

		// Read-only access
		float operator[](int idx) const {
			assert(idx >= 0 && idx < 3); // Ensure index is valid
			return (&x)[idx];
		}

		// Read/Write access
		float& operator[](int idx) {
			assert(idx >= 0 && idx < 3); // Ensure index is valid
			return (&x)[idx];
		}

		// --- Existing Methods ---

		// Note: This method causes a memory leak if the returned pointer is not deleted.
		float* toFloat3() {
			float* array = new float[3];
			array[0] = x;
			array[1] = y;
			array[2] = z;
			return array;
		}

		float dot(const Vec3& other) const {
			return x * other.x + y * other.y + z * other.z;
		}

		float length() const {
			return std::sqrt(x * x + y * y + z * z);
		}

		Vec3 normalized() const {
			float len = length();
			if (len > 1e-6f) {
				return { x / len, y / len, z / len };
			}
			return *this; // Avoid division by zero
		}

		Vec3 cross(const Vec3& other) const {
			return { y * other.z - z * other.y,
					z * other.x - x * other.z,
					x * other.y - y * other.x };
		}
	};

	// --- Non-Member Binary Operators (creates a new vector) ---

	// Vector addition
	inline Vec3 operator+(Vec3 lhs, const Vec3& rhs) {
		return lhs += rhs;
	}

	// Vector subtraction
	inline Vec3 operator-(Vec3 lhs, const Vec3& rhs) {
		return lhs -= rhs;
	}

	// Scalar multiplication
	inline Vec3 operator*(Vec3 lhs, float scalar) {
		return lhs *= scalar;
	}

	// Commutative scalar multiplication (allows float * Vec3)
	inline Vec3 operator*(float scalar, Vec3 rhs) {
		return rhs *= scalar;
	}

	// Scalar division
	inline Vec3 operator/(Vec3 lhs, float scalar) {
		return lhs /= scalar;
	}

	// Component-wise multiplication (Hadamard product)
	inline Vec3 operator*(Vec3 lhs, const Vec3& rhs) {
		lhs.x *= rhs.x;
		lhs.y *= rhs.y;
		lhs.z *= rhs.z;
		return lhs;
	}


	// --- Comparison Operators ---

	// Equality check
	inline bool operator==(const Vec3& lhs, const Vec3& rhs) {
		// Note: For robust floating point comparison, an epsilon check is often better.
		return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
	}

	// Inequality check
	inline bool operator!=(const Vec3& lhs, const Vec3& rhs) {
		return !(lhs == rhs);
	}

	// --- Stream Operator ---

	// For printing to streams like std::cout
	inline std::ostream& operator<<(std::ostream& os, const Vec3& v) {
		os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
		return os;
	}

struct Vertex
{
	// float* coords, * normals; //3d coordinates etc
	Vec3 coords,  normals; //3d coordinates etc
	int idx; //who am i; verts[idx]

	std::vector< int > vertList; //adj vertices;
	std::vector< int > triList; 
	std::vector< int > edgeList; 

	SbVec3f color;
	
//	Vertex(int i, float* c) : idx(i), coords(Vec3(c[0] ,c[1], c[2])) {};
	Vertex(int i, float x, float y, float z) : idx(i), coords(Vec3(x,y,z)) {};

};

struct Edge
{
	int idx; //edges[idx]
	int v1i, v2i; //endpnts
//	float length;   // USING EUCLIDIAN DISTANCE CALCULATION WHEN NECESSARY
  
	Edge(int id, int v1, int v2) : idx(id), v1i(v1), v2i(v2) {};

	// Overload the equality operator
	bool operator==(const Edge& other) const {
		// Consider edges (a, b) and (b, a) as the same
		return (v1i == other.v1i && v2i == other.v2i) ||
			(v1i == other.v2i && v2i == other.v1i);
	}
};

struct Triangle
{
	int idx; //tris[idx]
	int v1i, v2i, v3i;
  
	Triangle(int id, int v1, int v2, int v3) : idx(id), v1i(v1), v2i(v2), v3i(v3) {};
};

class Mesh
{
public:
	std::vector< Vertex* > vertices;
	std::vector< Triangle* > tris;
	std::vector< Edge* > edges;
	char* m_name;
	// painting
	std::vector< std::vector< Edge* > > highlightedEdges;
	std::vector< int > samples;
	std::vector< int > symm1;
	std::vector< int > symm2;
	std::vector< int > symmetryAxisVertices; // <-- ADDED THIS LINE

	Mesh();
	Mesh(char* name);
	~Mesh();
	// Copy constructor
	Mesh(Mesh& other);

	// Utility
	float getLength(int v1, int v2);
	float getTotalArea();
	void createCube(float side);
	void saveAsFile();
	void addTriangle(int v1, int v2, int v3);
	void normalizeMesh();
	void addEdge(int v1, int v2);
	void addVertex(float x, float y, float z);
	bool makeVertsNeighbor(int v1i, int v2i);

private:

};


} // namespace mesh
