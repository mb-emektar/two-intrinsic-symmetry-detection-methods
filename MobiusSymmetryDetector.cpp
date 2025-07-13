#include "MobiusSymmetryDetector.h"
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

MobiusSymmetryDetector::MobiusSymmetryDetector(mesh::Mesh* mesh, Dijkstra* a_dijkstra)
    : m_mesh(mesh), m_dijkstra(a_dijkstra) {
    if (!m_mesh) {
        throw std::invalid_argument("Mesh pointer cannot be null.");
    }
    if (!m_dijkstra) {
        throw std::invalid_argument("Dijkstra solver pointer cannot be null.");
    }
}

void MobiusSymmetryDetector::detectSymmetry() {
    if (m_mesh->vertices.empty()) {
        std::cerr << "Error: Mesh has no vertices." << std::endl;
        return;
    }
    if (m_mesh->samples.empty()) {
        // FIX: Added a fallback to use all vertices as samples if none are provided.
        std::cout << "Warning: No sample points provided in the mesh. Using all vertices as samples." << std::endl;
        m_mesh->samples.resize(m_mesh->vertices.size());
        for (int i = 0; i < m_mesh->vertices.size(); ++i) m_mesh->samples[i] = i;
    }

    // Step 1: Pre-computation of AGD.
    // Geodesic distances are calculated in runtime by the Dijkstra class.
    std::cout << "Calculating Average Geodesic Distance (AGD)..." << std::endl;
    calculateAGD();

    // Step 2: Flatten the mesh into the complex plane.
    std::cout << "Flattening mesh to complex plane with Conformal Map..." << std::endl;
    flattenMesh();
    //flattenMeshPCA();

    // Step 3: Search for the best anti-Möbius transformation
    std::cout << "Searching for best anti-Möbius transformation..." << std::endl;
    MobiusTransform best_transform = findBestAntiMobiusTransform(m_mesh->samples);

    // Step 4: Use the best transform to find and store the final correspondences.
    std::cout << "Extracting final correspondences..." << std::endl;
    extractCorrespondences(best_transform);

    std::cout << "Symmetry detection complete." << std::endl;
}


void MobiusSymmetryDetector::calculateAGD() {
    int numOfVertices = m_mesh->vertices.size();
    if (numOfVertices == 0) return;

    // --- IMPROVED AREA-WEIGHTED AGD CALCULATION ---

    // 1. Pre-calculate the area associated with each vertex.
    // The area of a vertex is 1/3 of the area of all incident triangles.
    // This correctly discretizes the surface integral in the AGD definition.
    std::vector<double> vertex_areas(numOfVertices, 0.0);
    for (const auto& tri : m_mesh->tris) {
        const mesh::Vec3& p1 = m_mesh->vertices[tri->v1i]->coords;
        const mesh::Vec3& p2 = m_mesh->vertices[tri->v2i]->coords;
        const mesh::Vec3& p3 = m_mesh->vertices[tri->v3i]->coords;

        // Area of a triangle = 0.5 * |(p2-p1) x (p3-p1)|
        float triangleArea = 0.5f * (p2 - p1).cross(p3 - p1).length();

        if (tri->v1i < numOfVertices) vertex_areas[tri->v1i] += triangleArea / 3.0;
        if (tri->v2i < numOfVertices) vertex_areas[tri->v2i] += triangleArea / 3.0;
        if (tri->v3i < numOfVertices) vertex_areas[tri->v3i] += triangleArea / 3.0;
    }

    m_agd.resize(numOfVertices);
    std::atomic<int> progressCounter(0);

    // Parallelize the calculation since each vertex's AGD is independent.
#pragma omp parallel for
    for (int i = 0; i < numOfVertices; ++i) {
        double agdSum = 0.0;

        // Get the geodesic distances from vertex 'i' to all other vertices.
        std::vector<float> distances = m_dijkstra->getDistanceVector(i);

        // 2. Calculate AGD as a discrete approximation of the integral:
        // Φ_agd(i) = ∫ d(i,j) d(area_j) ≈ Σ_j (distance(i,j) * area(j)) 
        for (int j = 0; j < numOfVertices; ++j) {
            agdSum += static_cast<double>(distances[j]) * vertex_areas[j];
        }
        m_agd[i] = agdSum;

        int completed = ++progressCounter;
        if (completed % (numOfVertices / 100 + 1) == 0 || completed == numOfVertices) {
#pragma omp critical
            {
                std::cout << "AGD Calculation Progress: " << (completed * 100) / numOfVertices << "%\r" << std::flush;
            }
        }
    }
    std::cout << std::endl;
}
void MobiusSymmetryDetector::flattenMeshPCA() {
    // PRINCIPAL COMPONENT ANALYSIS (PCA)

    int numOfVertices = m_mesh->vertices.size();
    if (numOfVertices < 3) {
        m_flattenedVertices.resize(numOfVertices);
        for (int i = 0; i < numOfVertices; ++i) {
            m_flattenedVertices[i] = std::complex<float>(m_mesh->vertices[i]->coords.x, m_mesh->vertices[i]->coords.y);
        }
        return;
    }

    // 1. Load vertex data into an Eigen matrix
    Eigen::MatrixXf points(numOfVertices, 3);
    for (int i = 0; i < numOfVertices; ++i) {
        points.row(i) << m_mesh->vertices[i]->coords.x, m_mesh->vertices[i]->coords.y, m_mesh->vertices[i]->coords.z;
    }

    // 2. Compute the mean of the vertices
    Eigen::RowVector3f centroid = points.colwise().mean();

    // 3. Center the points by subtracting the centroid
    points.rowwise() -= centroid;

    // 4. Compute the covariance matrix of the centered points
    Eigen::Matrix3f cov = (points.adjoint() * points) / float(numOfVertices - 1);

    // 5. Compute the eigenvectors of the covariance matrix.
    // SelfAdjointEigenSolver is efficient for symmetric matrices and sorts results by eigenvalue.
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigenSolver(cov);
    Eigen::Vector3f axis1 = eigenSolver.eigenvectors().col(2); // largest eigenvalue
    Eigen::Vector3f axis2 = eigenSolver.eigenvectors().col(1); // second largest eigenvalue

    // 6. Project the centered points onto the new 2D plane
    m_flattenedVertices.resize(numOfVertices);
    for (int i = 0; i < numOfVertices; ++i) {
        float u = points.row(i).dot(axis1);
        float v = points.row(i).dot(axis2);
        m_flattenedVertices[i] = std::complex<float>(u, v);
    }
}
void MobiusSymmetryDetector::flattenMesh() {
    // This is a simplier method instead of  a conformal parameterization.
    // The paper uses mid-edge uniformization .
    // I used a simple projection onto the XY plane, which is NOT conformal
    m_flattenedVertices.resize(m_mesh->vertices.size());
    for (size_t i = 0; i < m_mesh->vertices.size(); ++i) {
        m_flattenedVertices[i] = std::complex<float>(m_mesh->vertices[i]->coords.x, m_mesh->vertices[i]->coords.y);
    }
}
MobiusSymmetryDetector::MobiusTransform MobiusSymmetryDetector::findBestAntiMobiusTransform(const std::vector<int>& sampleIndices) {
    MobiusTransform bestTransform = MobiusTransform::Zero();
    float bestScore = -1;
    const float totalMeshArea = m_mesh->getTotalArea();
    if (totalMeshArea <= 0) return bestTransform;

    // A heuristic threshold to prune point pairs that are too close together.
    const float distPruneThreshold = std::sqrt(totalMeshArea / (32 * M_PI));
    std::atomic<int> progressCounter(0);
    const int numSamples = sampleIndices.size();

    // This is the main search loop, parallelized over the first point p1.
    // The complexity is O(n^3) on the number of samples n
#pragma omp parallel for
    for (int i = 0; i < numSamples; ++i) {
        int p1_idx = sampleIndices[i];

        std::vector<float> p1_distances = m_dijkstra->getDistanceVector(p1_idx);

        for (int j = i + 1; j < numSamples; ++j) {
            int p2_idx = sampleIndices[j];

            // 1: Skip pairs with very different AGD values.
            if (std::min(m_agd[p1_idx], m_agd[p2_idx]) / std::max(m_agd[p1_idx], m_agd[p2_idx]) < 0.9f) {
                continue;
            }
            // 2: Skip pairs that are too close on the mesh surface.
            if (p1_distances[p2_idx] < distPruneThreshold) {
                continue;
            }

            for (int k = 0; k < numSamples; ++k) {
                if (k == i || k == j) continue; // The three points must be distinct.
                int p_fixed_idx = sampleIndices[k];

                // Anti-Möbius transform is defined by 3 point correspondences.
                // Assuming (p1, p2) is a symmetric pair and p_fixed is on the symmetry axis.
                // So, m(conj(p1)) -> p2, m(conj(p2)) -> p1, and m(conj(p_fixed)) -> p_fixed.
                std::vector<std::pair<std::complex<float>, std::complex<float>>> constraints = {
                    { std::conj(m_flattenedVertices[p1_idx]), m_flattenedVertices[p2_idx] },
                    { std::conj(m_flattenedVertices[p2_idx]), m_flattenedVertices[p1_idx] },
                    { std::conj(m_flattenedVertices[p_fixed_idx]), m_flattenedVertices[p_fixed_idx] }
                };

                MobiusTransform candidateTransform = solveForMobius(constraints);
                if (candidateTransform.norm() < 1e-6) continue; // Invalid transform.

                float currsScore = measureAlignmentScore(candidateTransform);

                // Use a critical section to prevent a race condition when updating the best result.
#pragma omp critical
                {
                    if (currsScore > bestScore) {
                        bestScore = currsScore;
                        bestTransform = candidateTransform;
                    }
                }
            }
        }

        int completed = ++progressCounter;
        if (completed % (numSamples / 100 + 1) == 0 || completed == numSamples) {
#pragma omp critical
            {
                std::cout << "Search Progress: " << (completed * 100) / numSamples << "%\r" << std::flush;
            }
        }
    }
    std::cout << std::endl;
    return bestTransform;
}

bool MobiusSymmetryDetector::isPointInTriangle(const std::complex<float>& pt, const std::complex<float>& v1, const std::complex<float>& v2, const std::complex<float>& v3) const {
    auto sign = [](const std::complex<float>& p1, const std::complex<float>& p2, const std::complex<float>& p3) {
        return (p1.real() - p3.real()) * (p2.imag() - p3.imag()) - (p2.real() - p3.real()) * (p1.imag() - p3.imag());
        };

    float d1, d2, d3;
    bool has_neg, has_pos;

    d1 = sign(pt, v1, v2);
    d2 = sign(pt, v2, v3);
    d3 = sign(pt, v3, v1);

    has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
    has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

    return !(has_neg && has_pos);
}

// FIX: New helper function to calculate barycentric coordinates of a point within a 2D triangle.
mesh::Vec3 MobiusSymmetryDetector::barycentricCoordinates(const std::complex<float>& pt, const std::complex<float>& v1, const std::complex<float>& v2, const std::complex<float>& v3) const {
    float den = (v2.imag() - v3.imag()) * (v1.real() - v3.real()) + (v3.real() - v2.real()) * (v1.imag() - v3.imag());
    if (std::abs(den) < 1e-9) return mesh::Vec3(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0); // Degenerate triangle

    float w1 = ((v2.imag() - v3.imag()) * (pt.real() - v3.real()) + (v3.real() - v2.real()) * (pt.imag() - v3.imag())) / den;
    float w2 = ((v3.imag() - v1.imag()) * (pt.real() - v3.real()) + (v1.real() - v3.real()) * (pt.imag() - v3.imag())) / den;
    float w3 = 1.0f - w1 - w2;
    return mesh::Vec3(w1, w2, w3);
}

// FIX: This is the new, core function that correctly maps a 2D point back to the 3D surface.
mesh::Vec3 MobiusSymmetryDetector::map2DPointTo3DSurface(const std::complex<float>& point2D, int& found_tri_idx) const {
    // This is a simple but slow O(num_tris) implementation. For performance on large meshes,
    // a spatial data structure (like a quadtree or grid) on the flattened vertices would be needed.
    for (int i = 0; i < m_mesh->tris.size(); ++i) {
        const auto& tri = m_mesh->tris[i];
        const auto& z1 = m_flattenedVertices[tri->v1i];
        const auto& z2 = m_flattenedVertices[tri->v2i];
        const auto& z3 = m_flattenedVertices[tri->v3i];

        if (isPointInTriangle(point2D, z1, z2, z3)) {
            const auto& v1 = m_mesh->vertices[tri->v1i]->coords;
            const auto& v2 = m_mesh->vertices[tri->v2i]->coords;
            const auto& v3 = m_mesh->vertices[tri->v3i]->coords;

            mesh::Vec3 bary = barycentricCoordinates(point2D, z1, z2, z3);
            found_tri_idx = i;
            return v1 * bary.x + v2 * bary.y + v3 * bary.z;
        }
    }

    // Fallback: if point is not in any triangle, find the nearest vertex in the 2D plane
    // and return its 3D coordinates. This can happen due to numerical inaccuracies.
    float minDist_sq = std::numeric_limits<float>::max();
    int closestVert_idx = -1;
    for (size_t i = 0; i < m_flattenedVertices.size(); ++i) {
        float dist_sq = std::norm(m_flattenedVertices[i] - point2D);
        if (dist_sq < minDist_sq) {
            minDist_sq = dist_sq;
            closestVert_idx = i;
        }
    }
    found_tri_idx = -1; // Indicate no triangle was found
    return m_mesh->vertices[closestVert_idx]->coords;
}

float MobiusSymmetryDetector::measureAlignmentScore(const MobiusTransform& transform) {
    const std::vector<int>& sampleIndices = m_mesh->samples;
    if (sampleIndices.empty()) return 0.0f;

    const int numOfSamples = sampleIndices.size();
    std::vector<mesh::Vec3> original_3d_points(numOfSamples);
    std::vector<mesh::Vec3> transformed_3d_points(numOfSamples);

    // Step 1: Create the two sets of 3D points: P (originals) and P' (transformed)
#pragma omp parallel for
    for (int i = 0; i < numOfSamples; ++i) {
        int vertex_idx = sampleIndices[i];
        original_3d_points[i] = m_mesh->vertices[vertex_idx]->coords;

        std::complex<float> flat_pt = m_flattenedVertices[vertex_idx];
        std::complex<float> transformedFlat_pt = applyAntiMobius(transform, flat_pt);

        int tri_idx = -1; // Dummy variable
        transformed_3d_points[i] = map2DPointTo3DSurface(transformedFlat_pt, tri_idx);
    }

    std::vector<int> forwardMap(numOfSamples, -1);
    std::vector<int> backwardMap(numOfSamples, -1);

    // Step 2: For each point in P, find its closest neighbor in P' (Forward Pass)
    // This is the most expensive part.
#pragma omp parallel for
    for (int i = 0; i < numOfSamples; ++i) {
        float min_dist_sq = std::numeric_limits<float>::max();
        int best_j = -1;
        for (int j = 0; j < numOfSamples; ++j) {
            float dist_sq = (original_3d_points[i] - transformed_3d_points[j]).length();
            if (dist_sq < min_dist_sq) {
                min_dist_sq = dist_sq;
                best_j = j;
            }
        }
        forwardMap[i] = best_j;
    }

    // Step 3: For each point in P', find its closest neighbor in P (Backward Pass)
#pragma omp parallel for
    for (int j = 0; j < numOfSamples; ++j) {
        float min_dist_sq = std::numeric_limits<float>::max();
        int best_i = -1;
        for (int i = 0; i < numOfSamples; ++i) {
            float dist_sq = (transformed_3d_points[j] - original_3d_points[i]).length();
            if (dist_sq < min_dist_sq) {
                min_dist_sq = dist_sq;
                best_i = i;
            }
        }
        backwardMap[j] = best_i;
    }

    // Step 4: Count the mutually closest pairs
    int mutuallyClosestCount = 0;
    for (int i = 0; i < numOfSamples; ++i) {
        if (backwardMap[forwardMap[i]] == i) {
            mutuallyClosestCount++;
        }
    }

    // The score is the fraction of aligned points.
    return static_cast<float>(mutuallyClosestCount) / numOfSamples;
}


void MobiusSymmetryDetector::extractCorrespondences(const MobiusTransform& transform) {
    if (transform.norm() < 1e-6) {
        std::cerr << "Warning: Invalid (zero) transform provided. Cannot extract correspondences." << std::endl;
        return;
    }

    m_mesh->symm1.clear();
    m_mesh->symm2.clear();
    m_mesh->symmetryAxisVertices.clear();

    int numOfVertices = m_mesh->vertices.size();
    if (numOfVertices == 0) {
        return;
    }

    std::vector<int> correspondenceMap(numOfVertices, -1);

    // FIX: This section is also entirely rewritten to use the corrected 3D mapping.
#pragma omp parallel for
    for (int i = 0; i < numOfVertices; ++i) {
        std::complex<float> flat_pt = m_flattenedVertices[i];
        std::complex<float> transformed_flat_pt = applyAntiMobius(transform, flat_pt);

        // 1. Map the transformed 2D point back to a 3D point on the surface.
        int tri_idx = -1;
        mesh::Vec3 transformed_3d_pt = map2DPointTo3DSurface(transformed_flat_pt, tri_idx);

        // 2. Find the closest vertex on the entire mesh to this new 3D point.
        // This search is O(N), so this whole loop is O(N^2). Can be optimized with a 3D k-d tree.
        int closestVert_idx = -1;
        float minDist_sq = std::numeric_limits<float>::max();
        for (int j = 0; j < numOfVertices; ++j) {
            mesh::Vec3 P = m_mesh->vertices[j]->coords;
            float dist_sq = (P - transformed_3d_pt).length();
            if (dist_sq < minDist_sq) {
                minDist_sq = dist_sq;
                closestVert_idx = j;
            }
        }
        correspondenceMap[i] = closestVert_idx;
    }

    // Step 4: Identify symmetric pairs and symmetry axis vertices from the map.
    // This part is kept serial to avoid race conditions when modifying vectors.
    std::vector<bool> visited(numOfVertices, false);
    const float stationary_threshold = std::sqrt(m_mesh->getTotalArea() / (M_PI * 10000));
    int step = numOfVertices / 500;
    if (step == 0) step = 1;

    for (int i = 0; i < numOfVertices; ++i) {
        if (visited[i]) continue;

        int j = correspondenceMap[i];
        if (j == -1 || j >= numOfVertices) continue;

        // Case 1: Point is on the symmetry axis if it maps to itself.
        // We use the geodesic distance for a more robust check.
        float geodesic_dist = m_dijkstra->getDistanceVector(i)[j];
        if (geodesic_dist < stationary_threshold) {
            m_mesh->symmetryAxisVertices.push_back(i);
            visited[i] = true;
        }
        // Case 2: A pair is symmetric if they map to each other (are mutual correspondences).
        else if (correspondenceMap[j] == i) {
            if (i % step == 0) { // Only store a subset for visualization
                m_mesh->symm1.push_back(i);
                m_mesh->symm2.push_back(j);
            }
            visited[i] = true;
            visited[j] = true;
        }
    }
}

MobiusSymmetryDetector::MobiusTransform MobiusSymmetryDetector::solveForMobius(const std::vector<std::pair<std::complex<float>, std::complex<float>>>& pairs) {
    int k = pairs.size();
    if (k < 3) return MobiusTransform::Zero();

    Eigen::MatrixXcd A(k, 4);
    for (int i = 0; i < k; ++i) {
        std::complex<float> z = pairs[i].first;
        std::complex<float> w = pairs[i].second;
        // From az + b - cwz - dw = 0 
        A(i, 0) = z;
        A(i, 1) = 1;
        A(i, 2) = -w * z;
        A(i, 3) = -w;
    }

    // The solution to the homogeneous system Ax = 0 is the right singular vector
    // corresponding to the smallest singular value. This is the last column of V.
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(A, Eigen::ComputeFullV);
    return svd.matrixV().col(3);
}

std::complex<float> MobiusSymmetryDetector::applyMobius(const MobiusTransform& m, std::complex<float> z) {
    // The MobiusTransform uses double-precision complex numbers.
    // Cast them to float precision to match the input vertex data.
    std::complex<float> a(m(0).real(), m(0).imag());
    std::complex<float> b(m(1).real(), m(1).imag());
    std::complex<float> c(m(2).real(), m(2).imag());
    std::complex<float> d(m(3).real(), m(3).imag());

    if (std::abs(c * z + d) < 1e-9f) {
        return { std::numeric_limits<float>::infinity(), 0 };
    }
    return (a * z + b) / (c * z + d);
}

std::complex<float> MobiusSymmetryDetector::applyAntiMobius(const MobiusTransform& m, std::complex<float> z) {
    // An anti-Möbius transform is  actually a Möbius transform applied to the conj of z.
    return applyMobius(m, std::conj(z));
}