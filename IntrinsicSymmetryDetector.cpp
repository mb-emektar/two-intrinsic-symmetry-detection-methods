// IntrinsicSymmetryDetector.cpp
#include "IntrinsicSymmetryDetector.h"

namespace symm {

    IntrinsicSymmetryDetector::IntrinsicSymmetryDetector(Mesh* mesh)
        : m_mesh(mesh), m_kdTree(nullptr), m_annDataPoints(nullptr)
    {
        numOfVertices = static_cast<int>(m_mesh->vertices.size());
        mass_ = Eigen::VectorXd::Zero(numOfVertices);
    }

    IntrinsicSymmetryDetector::~IntrinsicSymmetryDetector() {
        cleanupKdTree();
        annClose();
    }

    void IntrinsicSymmetryDetector::cleanupKdTree() {
        if (m_kdTree) {
            delete m_kdTree;
            m_kdTree = nullptr;
        }
        if (m_annDataPoints) {
            annDeallocPts(m_annDataPoints);
            m_annDataPoints = nullptr;
        }
    }

    // Builds the kd-tree from the computed signatures for fast correspondence search.
    void IntrinsicSymmetryDetector::buildKdTree() {
        cleanupKdTree(); // Ensure no old tree exists.

        if (m_signatures.rows() == 0 || m_signatures.cols() == 0) {
            std::cerr << "Cannot build kd-tree: Signatures are empty." << std::endl;
            return;
        }

        int numPts = m_signatures.rows();
        int dim = m_signatures.cols();

        m_annDataPoints = annAllocPts(numPts, dim);
        for (int i = 0; i < numPts; ++i) {
            for (int j = 0; j < dim; ++j) {
                m_annDataPoints[i][j] = m_signatures(i, j);
            }
        }

        m_kdTree = new ANNkd_tree(m_annDataPoints, numPts, dim);
        std::cout << "Built kd-tree with " << numPts << " points of dimension " << dim << std::endl;
    }

    // Calculates the area of the triangle given its vertex indices.
    double IntrinsicSymmetryDetector::computeTriangleArea(int idx0, int idx1, int idx2) const {
        Vec3 v0 = m_mesh->vertices[idx0]->coords;
        Vec3 v1 = m_mesh->vertices[idx1]->coords;
        Vec3 v2 = m_mesh->vertices[idx2]->coords;
        Vec3 edge1 = v1 - v0;
        Vec3 edge2 = v2 - v0;
        return 0.5 * edge1.cross(edge2).length();
    }

    // Calculates the cotangent of the angle at vertex v0 in the triangle (v0, v1, v2).
    double IntrinsicSymmetryDetector::computeCotangent(const Vec3& v0, const Vec3& v1, const Vec3& v2) const {
        Vec3 a = v1 - v0;
        Vec3 b = v2 - v0;
        double dot_ab = a.dot(b);
        double crossLen = a.cross(b).length();

        // Avoid division by zero for degenerate triangles.
        if (crossLen < 1e-9) return 0;

        return dot_ab / crossLen;
    }

    // Calculates the cotangent Laplacian matrix L_ .
    void IntrinsicSymmetryDetector::buildLaplacianAndMass() {
        // 1. Calculates per-vertex mass. Each vertex gets 1/3 of the area of each adjacent triangle.
        for (auto tri : m_mesh->tris) {
            int i = tri->v1i;
            int j = tri->v2i;
            int k = tri->v3i;
            double area = computeTriangleArea(i, j, k);
            mass_(i) += area / 3;
            mass_(j) += area / 3;
            mass_(k) += area / 3;
        }

        // 2. Calculates cotangent weights for each edge.
        typedef Eigen::Triplet<double> T;
        std::vector<T> triplets;
        std::vector<std::unordered_map<int, double>> weights(numOfVertices);

        for (auto tri : m_mesh->tris) {
            int i = tri->v1i;
            int j = tri->v2i;
            int k = tri->v3i;
            Vec3 vi = m_mesh->vertices[i]->coords;
            Vec3 vj = m_mesh->vertices[j]->coords;
            Vec3 vk = m_mesh->vertices[k]->coords;

            double cot_i = computeCotangent(vi, vj, vk);
            double cot_j = computeCotangent(vj, vk, vi);
            double cot_k = computeCotangent(vk, vi, vj);

            weights[j][k] += cot_i;
            weights[k][j] += cot_i;
            weights[k][i] += cot_j;
            weights[i][k] += cot_j;
            weights[i][j] += cot_k;
            weights[j][i] += cot_k;
        }

        // 3. Assemble the sparse Laplacian matrix L from the weights.
        // Off-diagonal L(i,j) is -0.5 * (cot(alpha_ij) + cot(beta_ij)).
        // Diagonal L(i,i) is the sum of weights of all incident edges.
        std::vector<double> diag(numOfVertices, 0);
        for (int i = 0; i < numOfVertices; i++) {
            for (auto& entry : weights[i]) {
                int j = entry.first;
                double weight = 0.5 * entry.second;
                triplets.push_back(T(i, j, -weight));
                diag[i] += weight; 
            }
            triplets.push_back(T(i, i, diag[i]));
        }

        L_.resize(numOfVertices, numOfVertices);
        L_.setFromTriplets(triplets.begin(), triplets.end());

        // 4. Build the normalized Laplacian: A_norm = D^{-1/2} * L * D^{-1/2}, where D is the mass matrix.
        std::vector<T> tripletsNorm;
        for (int k = 0; k < L_.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(L_, k); it; ++it) {
                int i = it.row();
                int j = it.col();
                double value = it.value();
                double m_i = mass_(i);
                double m_j = mass_(j);
                if (m_i > 1e-9 && m_j > 1e-9) {
                    tripletsNorm.push_back(T(i, j, value / (std::sqrt(m_i) * std::sqrt(m_j))));
                }
                else {
                    tripletsNorm.push_back(T(i, j, value));
                }
            }
        }
        m_A_norm.resize(numOfVertices, numOfVertices);
        m_A_norm.setFromTriplets(tripletsNorm.begin(), tripletsNorm.end());
    }

    // Normalizes eigenvectors to have a consistent sign.
    // It makes subsequent computations deterministic.
    // It inverts any eigenvector whose element of largest magnitude is negative.
    void IntrinsicSymmetryDetector::canonizeEigenvectorSigns() {
        if (m_signatures.rows() == 0) return;

        for (int j = 0; j < m_signatures.cols(); ++j) {
            Eigen::Index maxAbsRow_idx;
            m_signatures.col(j).cwiseAbs().maxCoeff(&maxAbsRow_idx);

            if (m_signatures(maxAbsRow_idx, j) < 0) {
                m_signatures.col(j) *= -1;
            }
        }
    }

    // Performs eigen-decomposition on the normalized Laplacian to find the basis for the intrinsic space.
    bool IntrinsicSymmetryDetector::computeEigenDecomposition(int numEigs) {
        int numWanted = numEigs + 1;
        
        int ncv = std::min(2 * numWanted + 1, numOfVertices);

        // Use Spectra's shift-and-invert solver to efficiently find eigenvalues with the smallest magnitude,
        // as these correspond to the low-frequency shape information we're interested in.
        Spectra::SparseSymShiftSolve<double> op(m_A_norm);
        double sigma = -1e-9; // A small negative shift helps with the numerical stability of the solver.
        Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<double>> eigs(op, numWanted, ncv, sigma);

        std::cout << "Using robust shift-and-invert eigensolver to find smallest eigenvalues..." << std::endl;
        eigs.init();
        
        // The default number of iteration was 1000. However it was not enough for large models.
        int nconv = eigs.compute(Spectra::SortRule::LargestMagn, 10000, 1e-9);
        
        if (eigs.info() != Spectra::CompInfo::Successful) {
            return false;
        }

        Eigen::VectorXd evals = eigs.eigenvalues();
        Eigen::MatrixXd evecs = eigs.eigenvectors();

        // The solver returns eigenvalues with largest magnitude relative to the shift.
        // Reverse them to get smallest absolute magnitude first.
        evals.reverseInPlace();
        evecs = evecs.rowwise().reverse().eval();

        // Create the signatures matrix. Skip the first trivial eigenvector (j=0).
        // The signature is the eigenvector scaled by 1/sqrt(eigenvalue).
        int d = numEigs;
        eigenvalues_.resize(d);
        m_signatures.resize(numOfVertices, d);
        if (numOfVertices == 0) return true;

        for (int j = 1; j < numWanted; j++) { 
            int signatureCol_idx = j - 1;   
            double lambda = evals(j);
            eigenvalues_(signatureCol_idx) = lambda;
            Eigen::VectorXd phi = evecs.col(j);

            double scalingFactor = (std::abs(lambda) < 1e-6) ? 0 : 1 / std::sqrt(std::abs(lambda));
            m_signatures.col(signatureCol_idx) = phi * scalingFactor;
        }

        // Make the signatures robust and consistent with post-processing.
        canonizeEigenvectorSigns();
        findOptimalEigenbasis();

        // Finally, build the kd-tree usin the processed signatures.
        buildKdTree();

        return true;
    }

    // Rotates the eigenbasis using SVD to make the signatures robust against pose changes.
    // This finds a canonical basis aligned with the principal components of the signature data.
    void IntrinsicSymmetryDetector::findOptimalEigenbasis() {
        if (m_signatures.rows() == 0) return;

        // 1. Compute the covariance-like matrix C = signatures^T * signatures.
        Eigen::MatrixXd C = m_signatures.transpose() * m_signatures;

        // 2. Perform SVD on C. The columns of U are the principal directions.
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(C, Eigen::ComputeThinU);
        Eigen::MatrixXd R = svd.matrixU(); // This is the rotation matrix.

        // 3. Apply the rotation to the original signatures.
        m_signatures = m_signatures * R;

        // 4. Re-canonize signs after rotation to ensure consistency.
        canonizeEigenvectorSigns();

        std::cout << "Applied SVD-based rotation to find an optimal eigenbasis." << std::endl;
    }

    // Finds the symmetric correspondence for every vertex
    std::vector<int> IntrinsicSymmetryDetector::computeSymmetryCorrespondences(const Eigen::VectorXi& S) const {
        if (!m_kdTree) {
            std::cerr << "Error: kd-tree not built. Call computeEigenDecomposition first." << std::endl;
            return std::vector<int>(numOfVertices, -1);
        }
        if (S.size() != m_signatures.cols()) {
            std::cerr << "Dimension mismatch: sign sequence size does not match signature dimension." << std::endl;
            return std::vector<int>();
        }

        std::vector<int> correspondences(numOfVertices, -1);
        int d = S.size();

        ANNpoint queryPoint = annAllocPt(d);
        ANNidx nnIdx;
        ANNdist distSq;

        for (int i = 0; i < numOfVertices; i++) {
            Eigen::VectorXd sPrime = m_signatures.row(i).transpose();
            for (int j = 0; j < d; j++) {
                sPrime(j) *= S(j);
            }

            for (int dim_idx = 0; dim_idx < d; ++dim_idx) {
                queryPoint[dim_idx] = sPrime(dim_idx);
            }

            m_kdTree->annkSearch(queryPoint, 1, &nnIdx, &distSq, 0);
            correspondences[i] = nnIdx;
        }

        annDeallocPt(queryPoint);
        return correspondences;
    }

    // Calculates the total squared Euclidean distance error for a given symmetry sign sequence.
    // This is the measure of how "good" the symmetry is. Smaller the error better the symmetry.
    double IntrinsicSymmetryDetector::computeErrorForSignSequence(const Eigen::VectorXi& S) const {
        if (!m_kdTree) {
            std::cerr << "Error: kd-tree not built." << std::endl;
            return std::numeric_limits<double>::infinity();
        }

        int d = S.size();
        if (d != m_signatures.cols()) {
            std::cerr << "Dimension mismatch in computeErrorForSignSequence." << std::endl;
            return std::numeric_limits<double>::infinity();
        }

        double totalError_E_S = 0;
        ANNpoint queryPoint = annAllocPt(d);
        ANNidx nnIdx;
        ANNdist distSq;

        for (int q_idx = 0; q_idx < numOfVertices; ++q_idx) {
            Eigen::VectorXd sPrime_q = m_signatures.row(q_idx).transpose();
            for (int j = 0; j < d; ++j) {
                sPrime_q(j) *= S(j);
            }

            for (int dim_idx = 0; dim_idx < d; ++dim_idx) {
                queryPoint[dim_idx] = sPrime_q(dim_idx);
            }

            m_kdTree->annkSearch(queryPoint, 1, &nnIdx, &distSq, 0);
            totalError_E_S += distSq; 
        }

        annDeallocPt(queryPoint);
        return totalError_E_S;
    }

    // A brute-force approach to find the best sign sequences by testing every possible combination.
    void IntrinsicSymmetryDetector::generateSignSequencesBruteForce(int numTopSequencesToKeep) {
        bestSignSequences_.clear();
        if (m_signatures.cols() == 0 || !m_kdTree) {
            std::cerr << "Cannot generate sequences: Signatures or kd-tree not ready." << std::endl;
            return;
        }

        int d = m_signatures.cols();
        int practical_d = std::min(d, 20);
        if (d > practical_d) {
            std::cout << "Warning: Brute-force search capped at d=" << practical_d << " (original d=" << d << ")" << std::endl;
        }

        long long numTotalSequences = 1LL << practical_d;
        std::cout << "Brute-force generating " << numTotalSequences << " sign sequences..." << std::endl;

        std::vector<SignSequenceWithError> allSequencesWithErrors;
        for (long long i = 0; i < numTotalSequences; ++i) {
            Eigen::VectorXi current_S = Eigen::VectorXi::Ones(d);
            for (int j = 0; j < practical_d; ++j) {
                if ((i >> j) & 1) {
                    current_S(j) = -1;
                }
            }
            double error = computeErrorForSignSequence(current_S);
            allSequencesWithErrors.push_back({ current_S, error });
        }

        std::cout << "Sorting sequences by error..." << std::endl;
        std::sort(allSequencesWithErrors.begin(), allSequencesWithErrors.end());

        int count_to_keep = std::min((int)allSequencesWithErrors.size(), numTopSequencesToKeep);
        for (int i = 0; i < count_to_keep; ++i) {
            bestSignSequences_.push_back(allSequencesWithErrors[i].sequence);
        }
    }

    // Generates candidate sign sequences using a multi-stage heuristic avoiding the brute-force search.
    void IntrinsicSymmetryDetector::generateSignSequences(int numTopSequencesToKeep) {
        bestSignSequences_.clear();
        if (m_signatures.cols() == 0 || numOfVertices == 0) {
            std::cerr << "Signatures have not been computed; cannot generate sequences." << std::endl;
            return;
        }
        int d = m_signatures.cols();

        std::cout << "Generating sign sequences using heuristic method..." << std::endl;

        // Stage 1: Detect eigenfunctions with low error when their sign is inverted individually.
        std::vector<std::pair<double, int>> potentiallyNegativeInfo = detectPotentiallyNegativeEigenfunctions();

        if (potentiallyNegativeInfo.empty()) {
            std::cout << "Heuristic did not find any potentially negative eigenfunctions." << std::endl;
        }

        // Stage 2: Group the candidate eigenfunctions.
        std::vector<std::vector<int>> coupledGroups = groupCoupledEigenfunctions(potentiallyNegativeInfo);

        // Stage 3: Form the final sign sequences from the groups.
        bestSignSequences_.push_back(Eigen::VectorXi::Ones(d));

        for (const auto& group : coupledGroups) {
            Eigen::VectorXi S = Eigen::VectorXi::Ones(d);
            for (int sigIdx : group) {
                if (sigIdx >= 0 && sigIdx < d) {
                    S(sigIdx) = -1;
                }
            }
            bestSignSequences_.push_back(S);
        }

        // Sort the sequences by error and keep only the best ones.
        if (numTopSequencesToKeep > 0 && bestSignSequences_.size() > numTopSequencesToKeep) {
            if (!m_kdTree) {
                std::cerr << "kd-tree not built, cannot rank heuristic sequences by error. Keeping all found." << std::endl;
            }
            else {
                std::vector<SignSequenceWithError> sequencesWithErrors;
                for (const auto& s_seq : bestSignSequences_) {
                    sequencesWithErrors.push_back({ s_seq, computeErrorForSignSequence(s_seq) });
                }
                std::sort(sequencesWithErrors.begin(), sequencesWithErrors.end());

                bestSignSequences_.clear();
                int count_to_keep = std::min((int)sequencesWithErrors.size(), numTopSequencesToKeep);
                for (int i = 0; i < count_to_keep; ++i) {
                    bestSignSequences_.push_back(sequencesWithErrors[i].sequence);
                }
            }
        }
        std::cout << "Finished generating sign sequences via heuristic. Found " << bestSignSequences_.size() << " sequences." << std::endl;
    }

    // build a temporary kd-tree from an arbitrary matrix of points.
    ANNkd_tree* IntrinsicSymmetryDetector::buildTemporaryKdTree(
        const Eigen::MatrixXd& points, ANNpointArray& annPointsStorage) const {

        if (points.rows() == 0) {
            annPointsStorage = nullptr;
            return nullptr;
        }
        int numPts = points.rows();
        int dim = points.cols();

        annPointsStorage = annAllocPts(numPts, dim);
        for (int i = 0; i < numPts; ++i) {
            for (int j = 0; j < dim; ++j) {
                annPointsStorage[i][j] = points(i, j);
            }
        }
        return new ANNkd_tree(annPointsStorage, numPts, dim);
    }

    // 1: Identify "potentially negative" eigenfunctions.
    // If the error from inverting its sign is small, an eigenfunction is negative 
    std::vector<std::pair<double, int>> IntrinsicSymmetryDetector::detectPotentiallyNegativeEigenfunctions() const {
        std::vector<std::pair<double, int>> potentialNegativesWithErrors;
        if (m_signatures.cols() == 0) return potentialNegativesWithErrors;

        int d = m_signatures.cols();
        Eigen::MatrixXd sPlusDataset(numOfVertices, d);
        ANNpointArray ann_s_plusPoints = nullptr;
        ANNkd_tree* tempKdTree = nullptr;
        ANNpoint queryPoint = annAllocPt(d);

        for (int k_idx = 0; k_idx < d; ++k_idx) { 
            // 1. Prepare the dataset where all signatures have positive components
            for (int q_vtx = 0; q_vtx < numOfVertices; ++q_vtx) {
                for (int comp = 0; comp < d; ++comp) {
                    sPlusDataset(q_vtx, comp) = (comp == k_idx) ? m_signatures(q_vtx, comp) : std::abs(m_signatures(q_vtx, comp));
                }
            }

            // 2. Build kd-tree on this modified set
            if (tempKdTree) delete tempKdTree;
            if (ann_s_plusPoints) annDeallocPts(ann_s_plusPoints);
            tempKdTree = buildTemporaryKdTree(sPlusDataset, ann_s_plusPoints);
            if (!tempKdTree) continue;

            // 3. Compute the error
            double E_k = 0;
            ANNidx nnIdx;
            ANNdist distSq;
            for (int pVtx = 0; pVtx < numOfVertices; ++pVtx) {
                for (int comp = 0; comp < d; ++comp) {
                    queryPoint[comp] = (comp == k_idx) ? -m_signatures(pVtx, comp) : std::abs(m_signatures(pVtx, comp));
                }
                tempKdTree->annkSearch(queryPoint, 1, &nnIdx, &distSq, 0);
                E_k += distSq;
            }
            potentialNegativesWithErrors.push_back({ E_k, k_idx });
        }

        if (tempKdTree) delete tempKdTree;
        if (ann_s_plusPoints) annDeallocPts(ann_s_plusPoints);
        if (queryPoint) annDeallocPt(queryPoint);

        // Sort by error and find the the smallest error.
        std::sort(potentialNegativesWithErrors.begin(), potentialNegativesWithErrors.end());

        int numToConsider = std::max(1, d / 2);
        if (potentialNegativesWithErrors.size() > numToConsider) {
            potentialNegativesWithErrors.resize(numToConsider);
        }
        return potentialNegativesWithErrors;
    }

    // Checks if two eigenfunctions are "coupled".
    // They are coupled if inverting one of them but not the other results with a high error.
    bool IntrinsicSymmetryDetector::areEigenfunctionsCoupled(
        int sigIdx1, int sigIdx2, double E_sigIdx1, double E_sigIdx2) const {

        if (m_signatures.cols() == 0) return false;
        int d = m_signatures.cols();

        Eigen::MatrixXd s_ij_plus_dataset(numOfVertices, d);
        ANNpointArray ann_s_ij_plus_points = nullptr;
        ANNkd_tree* tempKdTree = nullptr;
        ANNpoint queryPoint = annAllocPt(d);

        for (int q_vtx = 0; q_vtx < numOfVertices; ++q_vtx) {
            for (int comp = 0; comp < d; ++comp) {
                s_ij_plus_dataset(q_vtx, comp) = (comp == sigIdx1 || comp == sigIdx2) ? m_signatures(q_vtx, comp) : std::abs(m_signatures(q_vtx, comp));
            }
        }
        tempKdTree = buildTemporaryKdTree(s_ij_plus_dataset, ann_s_ij_plus_points);
        if (!tempKdTree) {
            if (queryPoint) annDeallocPt(queryPoint);
            return false;
        }

        double E_ij = 0, E_ji = 0;
        ANNidx nnIdx;
        ANNdist distSq;

        for (int p_vtx = 0; p_vtx < numOfVertices; ++p_vtx) {
            for (int comp = 0; comp < d; ++comp) queryPoint[comp] = s_ij_plus_dataset(p_vtx, comp);
            queryPoint[sigIdx2] *= -1;
            tempKdTree->annkSearch(queryPoint, 1, &nnIdx, &distSq, 0);
            E_ij += distSq;
        }

        for (int p_vtx = 0; p_vtx < numOfVertices; ++p_vtx) {
            for (int comp = 0; comp < d; ++comp) queryPoint[comp] = s_ij_plus_dataset(p_vtx, comp);
            queryPoint[sigIdx1] *= -1;
            tempKdTree->annkSearch(queryPoint, 1, &nnIdx, &distSq, 0);
            E_ji += distSq;
        }

        if (tempKdTree) delete tempKdTree;
        if (ann_s_ij_plus_points) annDeallocPts(ann_s_ij_plus_points);
        if (queryPoint) annDeallocPt(queryPoint);

        // Coupling condition: If flipping one but not the other does not
        // result in a small error, they are probably coupled.
        double zeroThreshold = 1e-5 * numOfVertices;
        if (E_ij < zeroThreshold || E_ji < zeroThreshold) {
            return false; 
        }

        // Coupled if the mixed-sign error is worse than the individual negative errors.
        bool coupled = (E_ij > std::min(E_sigIdx1, E_sigIdx2)) && (E_ji > std::min(E_sigIdx1, E_sigIdx2));
        return coupled;
    }

    // Groups the coupled eigenfunctions.
    // Each group represents a set of eigenfunctions that must invert their signs together.
    std::vector<std::vector<int>> IntrinsicSymmetryDetector::groupCoupledEigenfunctions(
        const std::vector<std::pair<double, int>>& potentiallyNegativeInfo) const {

        std::vector<std::vector<int>> groups;
        if (potentiallyNegativeInfo.empty()) return groups;

        std::vector<int> pniIndices;
        std::map<int, double> pniErrorsMap;
        for (const auto& p : potentiallyNegativeInfo) {
            pniIndices.push_back(p.second);
            pniErrorsMap[p.second] = p.first;
        }

        int numOfPni = pniIndices.size();
        std::vector<bool> visited(numOfPni, false);
        std::vector<std::vector<int>> adj(numOfPni);

        for (int i = 0; i < numOfPni; ++i) {
            for (int j = i + 1; j < numOfPni; ++j) {
                int sigIdx1 = pniIndices[i];
                int sigIdx2 = pniIndices[j];
                if (areEigenfunctionsCoupled(sigIdx1, sigIdx2, pniErrorsMap.at(sigIdx1), pniErrorsMap.at(sigIdx2))) {
                    adj[i].push_back(j);
                    adj[j].push_back(i);
                }
            }
        }

        // Finds the connected components in the graph
        for (int i = 0; i < numOfPni; ++i) {
            if (!visited[i]) {
                std::vector<int> currentgroup;
                std::queue<int> q;
                q.push(i);
                visited[i] = true;
                while (!q.empty()) {
                    int u_local = q.front(); q.pop();
                    currentgroup.push_back(pniIndices[u_local]);
                    for (int vLocal : adj[u_local]) {
                        if (!visited[vLocal]) {
                            visited[vLocal] = true;
                            q.push(vLocal);
                        }
                    }
                }
                groups.push_back(currentgroup);
            }
        }
        return groups;
    }

    // Finds the longest axis of symmetry from the set of fixed points.
    std::vector<int> IntrinsicSymmetryDetector::findLongestSymmetryAxis(const std::vector<int>& correspondences) const {
        if (!m_mesh || m_mesh->vertices.empty()) return {};

        // 1. Find all fixed points
        std::unordered_set<int> fixedPoints;
        for (int i = 0; i < correspondences.size(); ++i) {
            if (correspondences[i] == i) {
                fixedPoints.insert(i);
            }
        }
        if (fixedPoints.empty()) return {};

        // 2. Find all connected components of the fixed points with BFS.
        std::vector<std::vector<int>> components;
        std::unordered_set<int> visited;

        for (int startNode : fixedPoints) {
            if (visited.find(startNode) == visited.end()) {
                std::vector<int> currentComponent;
                std::queue<int> q;
                q.push(startNode);
                visited.insert(startNode);

                while (!q.empty()) {
                    int u = q.front(); q.pop();
                    currentComponent.push_back(u);

                    // A neighbor is in the same component if it's also a fixed point.
                    for (int v : m_mesh->vertices[u]->vertList) {
                        if (fixedPoints.count(v) && visited.find(v) == visited.end()) {
                            visited.insert(v);
                            q.push(v);
                        }
                    }
                }
                components.push_back(currentComponent);
            }
        }

        // 3. Return the largest component. It represents the longest continuous axis.
        if (components.empty()) return {};

        size_t longest_idx = 0;
        for (size_t i = 1; i < components.size(); ++i) {
            if (components[i].size() > components[longest_idx].size()) {
                longest_idx = i;
            }
        }
        return components[longest_idx];
    }

} // namespace symm