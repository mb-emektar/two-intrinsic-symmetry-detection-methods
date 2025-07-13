#pragma once

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues> 
#include <vector>
#include <complex>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <queue>
#include <atomic>
#include <omp.h>

#include "Mesh.h"
#include "Dijkstra.h"

/**
 * @class MobiusSymmetryDetector
 * @brief Detects reflective symmetries in a 3D mesh using Möbius transformation.
 *
 * This class implements the symmetry detection pipeline described in "Symmetry
 * in 3D Geometry: Extraction and Applications" by Mitra et al. The core idea is
 * flattening the mesh to the complex plane. Then, find an anti-Möbius transformation
 * that best maps the flattened shape to itself. Finally, use this transformation
 * to identify symmetric vertices on the original 3D mesh.
 *
 */
class MobiusSymmetryDetector {
public:
    /**
     * @brief Constructs the symmetry detector.
     * @param mesh: pointer to the mesh object
     * @param dijkstraSolver: pointer to a Dijkstra
     */
    explicit MobiusSymmetryDetector(mesh::Mesh* mesh, Dijkstra* a_dijkstra);

    /**
     * @brief Runs the main symmetry detection pipeline.
     *
     * Here are the steps of the function:
     * 1. Computes the Average Geodesic Distance (AGD) for vertices.
     * 2. Flattens the mesh to the 2D complex plane (a simplified uniformization).
     * 3. Searches for the optimal anti-Möbius transformation using sample points.
     * 4. Extracts vertex correspondences (symmetric pairs and stationary points)
     * based on the best transformation found.
     */
    void detectSymmetry();

private:
    // A Möbius transformation is defined by four complex coefficients m(z) = (az+b)/(cz+d).
    using MobiusTransform = Eigen::Vector4cd;

    /**
     * @brief Calculates the Average Geodesic Distance (AGD) for each vertex.
     * It is a measure of how "central" a vertex is
     */
    void calculateAGD();

    /**
     * @brief Flattens the 3D mesh onto the 2D complex plane.
     * It is not conformal mapping. It is the much simpler way to flatten a mesh.
	 * Since conformal mapping is really expensive to compute, I used a simple 
     * but effective method.
     */
    void flattenMesh();

    /**
    * @brief Flattens the 3D mesh onto the 2D complex plane.
    * It is PCA based method.
    */
    void flattenMeshPCA();

    /**
     * @brief Searches for the best anti-Möbius transformation and finds the one
     * maps the mesh to itself.
     * @param sampleIndices: vector of vertex indices to use for generating candidates.
     * @return The best anti-Möbius transformation
     */
    MobiusTransform findBestAntiMobiusTransform(const std::vector<int>& sampleIndices);

    /**
     * @brief Measures how well a given transformation aligns the surface with itself.
     * @param transform: The anti-Möbius transformation to score.
     * @return The alignment score, from 0 as no alignment to 1 as perfect alignment.
     */
    float measureAlignmentScore(const MobiusTransform& transform);

    /**
     * @brief Extracts final symmetric correspondences for all vertices.
     * @param transform : The best anti-Möbius transformation.
     */
    void extractCorrespondences(const MobiusTransform& transform);

    /**
     * @brief Solves for the coefficients of a Möbius transformation from point pairs.
     * @param pairs : A vector of corresponding point pairs (z, w).
     * @return The Möbius transformation coefficients .
     */
    MobiusTransform solveForMobius(const std::vector<std::pair<std::complex<float>, std::complex<float>>>& pairs);

    /**
     * @brief Conducts a Möbius transformation -> m(z) = (az+b)/(cz+d) to a point.
     * @param m The transformation coefficients [a, b, c, d].
     * @param z The complex number to transform.
     * @return The transformed complex number.
     */
    std::complex<float> applyMobius(const MobiusTransform& m, std::complex<float> z);

    /**
     * @brief Applies an anti-Möbius transformation m(conj(z)) to a point.
     * @param m The transformation coefficients.
     * @param z The complex number to transform.
     * @return The transformed complex number.
     */
    std::complex<float> applyAntiMobius(const MobiusTransform& m, std::complex<float> z);

    // --- Helpers for mapping 2D points back to the 3D surface ---

    /**
     * @brief Maps a 2D point from the flattened plane back to a 3D point on the mesh surface.
     * @param point2D The 2D point to map.
     * @param found_tri_idx A reference to store the index of the triangle containing the point.
     * @return The corresponding 3D point on the mesh surface.
     */
    mesh::Vec3 map2DPointTo3DSurface(const std::complex<float>& point2D, int& found_tri_idx) const;

    /**
     * @brief Calculates the barycentric coordinates of a point within a 2D triangle.
     * @return A Vec3 containing the barycentric coordinates (u, v, w).
     */
    mesh::Vec3 barycentricCoordinates(const std::complex<float>& pt, const std::complex<float>& v1, const std::complex<float>& v2, const std::complex<float>& v3) const;

    /**
     * @brief Determines if a 2D point is inside a 2D triangle.
     * @return True if the point is inside the triangle, false otherwise.
     */
    bool isPointInTriangle(const std::complex<float>& pt, const std::complex<float>& v1, const std::complex<float>& v2, const std::complex<float>& v3) const;

    Dijkstra* m_dijkstra;   
    mesh::Mesh* m_mesh;     

    std::vector<float> m_agd; // Stores the calculated AGD 
    std::vector<std::complex<float>> m_flattenedVertices; // 2D representation of the vertices
};