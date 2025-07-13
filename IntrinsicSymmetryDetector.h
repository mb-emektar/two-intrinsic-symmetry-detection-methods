// IntrinsicSymmetryDetector.h
#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>

#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>

#include "Mesh.h"

#include <vector>
#include <unordered_map>
#include <cmath>
#include <iostream>
#include <limits>
#include <algorithm> 
#include <queue>
#include <map> 
#include <unordered_set>

#include <ANN/ANN.h> 

using namespace Eigen;
using namespace Spectra;
using namespace mesh;

namespace symm {

    /**
     * @class IntrinsicSymmetryDetector
     * @brief Detects intrinsic symmetries of a mesh.
     *
     * This class implements the approach represented in "Intrinsic Symmetry Detection"
     * by Ovsjanikov et al.
     */
    class IntrinsicSymmetryDetector {
    public:
        /**
         * @brief Constructor
         * @param mesh A pointer to the Mesh object.
         */
        IntrinsicSymmetryDetector(Mesh* mesh);

        /**
         * @brief Destructor. Cleans up dynamically allocated resources, like the ANN kd-tree.
         */
        ~IntrinsicSymmetryDetector();

        // --- Core Pipeline Methods ---

        /**
         * @brief Builds the cotangent Laplacian matrix (L).
         */
        void buildLaplacianAndMass();

        /**
         * @brief Computes the eigen-decomposition of the normalized Laplacian.
         * The eigenvalues and eigenvectors form the basis for the intrinsic signatures.
         * @param numEigs The number of non-trivial eigenpairs to compute defining the dimension of the signature space.
         * @return True if the computation was successful, false otherwise.
         */
        bool computeEigenDecomposition(int numEigs);

        /**
         * @brief Generates candidate sign sequences that contain potential symmetries.
         * @param numTopSequencesToKeep The number of the best candidate sequences to store.
         */
        void generateSignSequences(int numTopSequencesToKeep = 5);

        /**
         * @brief Finds its corresponding symmetric vertex based on a given sign sequence for each vertex.
         * @param S The sign sequence which represents the symmetry transformation.
         * @return A vector showing each element at index `i` is the index of the vertex symmetric to vertex `i`.
         */
        std::vector<int> computeSymmetryCorrespondences(const Eigen::VectorXi& S) const;

        /**
         * @brief Finds the longest continuous path of vertices that are associated to the symmetry.
         * These paths represent the geometric axes of symmetry.
         * @param The map of symmetry correspondences per vertex.
         * @return A vector of vertex indices that represents the longest symmetry axis.
         */
        std::vector<int> findLongestSymmetryAxis(const std::vector<int>& correspondences) const;

        // --- Getters and OTher Utility Functions ---

        /**
         * @brief Gets the computed intrinsic signatures.
         * @return A constant reference to the signature matrix. Each row is the signature for a vertex.
         */
        const Eigen::MatrixXd& getSignatures() const { return m_signatures; }

        /**
         * @brief Gets the candidate sign sequences found with the `generateSignSequences` call.
         * @return A constant reference to the vector of sign sequences.
         */
        const std::vector<Eigen::VectorXi>& getSignSequences() const { return bestSignSequences_; }

        /**
         * @brief Copmutes the overall error for a particular sign sequence.
         * The error measures how well the given sequence represents a valid symmetry. The samller error is better.
         * @param S The sign sequence to evaluate.
         * @return The error value.
         */
        double computeErrorForSignSequence(const Eigen::VectorXi& S) const;

        /**
         * @brief A brute-force method to generate sign sequences by testing one by one.
         * Tests all 2^d possible sequences.
         * BEWARE: This method is computationally expensive and should only be used for small signature dimensions (d < 20).
         * @param numTopSequencesToKeep The number of best sequences to store which are ranked by their error.
         */
        void generateSignSequencesBruteForce(int numTopSequencesToKeep = 5);

    private:
        // --- Core Member Data ---
        Mesh* m_mesh;
        int numOfVertices;

        // --- Spectral Analysis Data ---
        Eigen::SparseMatrix<double> L_;         // Cotangent Laplacian matrix
        Eigen::SparseMatrix<double> m_A_norm;    // Normalized Laplacian matrix
        Eigen::VectorXd mass_;                  // Per-vertex mass/area vector
        Eigen::MatrixXd m_signatures;            // Intrinsic signatures (scaled eigenvectors)
        Eigen::VectorXd eigenvalues_;           // Eigenvalues from decomposition

        // --- Symmetry Data ---
        std::vector<Eigen::VectorXi> bestSignSequences_;

        // Helper struct to sort sign sequences by their associated error.
        struct SignSequenceWithError {
            Eigen::VectorXi sequence;
            double error;
            bool operator<(const SignSequenceWithError& other) const {
                return error < other.error;
            }
        };

        // --- ANN kd-tree for Nearest Neighbor Search ---
        ANNkd_tree* m_kdTree;        
        ANNpointArray m_annDataPoints;

        // --- Private Helper Methods ---

        double computeTriangleArea(int idx0, int idx1, int idx2) const;
        double computeCotangent(const Vec3& v0, const Vec3& v1, const Vec3& v2) const;

        void buildKdTree();
        void cleanupKdTree();

        void canonizeEigenvectorSigns();
        void findOptimalEigenbasis();

        // --- Heuristic Method for Sign Sequence Generation ---

        // 1: Identify eigenfunctions which are possibly be part of a symmetry.
        std::vector<std::pair<double, int>> detectPotentiallyNegativeEigenfunctions() const;

        // 2: Group the eigenfunctions that must invert their signs together.
        std::vector<std::vector<int>> groupCoupledEigenfunctions(
            const std::vector<std::pair<double, int>>& potentiallyNegativeInfo) const;

        // Helper for Stage 2: Checks if any are there any same eigenvectors.
        bool areEigenfunctionsCoupled(int sigIdx1, int sigIdx2,
            double E_sigIdx1, double E_sigIdx2) const;

        // build a temporary kd-tree for the heuristic steps.
        ANNkd_tree* buildTemporaryKdTree(const Eigen::MatrixXd& points, ANNpointArray& annPoints) const;
    };

} // namespace symm