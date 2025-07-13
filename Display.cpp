#include "Display.h"
#include <thread>
#include <chrono>   // For duration calculation
#include <fstream>  // For file output
#include <iomanip>  // For output formatting (std::fixed, std::setprecision)


Display::Display(char* argv, char* name)
{
    m_meshName = name;
    window = SoWin::init(argv);

    viewer = new SoWinExaminerViewer(window);

    root = new SoSeparator;
    root->ref();

    m_mesh = new Mesh(m_meshName);

    m_mesh->normalizeMesh(); // Normalize the mesh before processing
    m_painter = new Painter();
}

Display::~Display()
{
}

void Display::displayMesh()
{
    root->addChild(m_painter->getShapeSep(m_mesh));

    viewer->setSize(SbVec2s(1080, 1080));
    viewer->setSceneGraph(root);
    viewer->show();
    SoWin::show(window);
    SoWin::mainLoop();
    delete viewer;
    root->unref();
    SoWin::exitMainLoop();
}

void Display::createMembers(int methodSelection)
{
    m_dijsktra = new Dijkstra(m_mesh);
    if(methodSelection==2)
        m_dijsktra->precomputeAndCacheAllDistances();
}

void Display::calculateSampling()
{
    std::cout << "Please enter the number of sample points: " << std::endl;
    std::cin >> numOfSamplePts;
    if (numOfSamplePts < 0)
    {
        numOfSamplePts = 0;
    }

    std::cout << "Calculating Farthest Point Sampling..." << std::endl;
    m_sampling = new Sampling(m_mesh, m_dijsktra, numOfSamplePts); //manageri ve mesh i ayri ayri ver
    m_sampling->calculateSamples();
    std::cout << "Sampling Completed..." << std::endl;
}
// Integration into the Display class
void Display::calculateMobiusSymmetry() {
    auto start = std::chrono::high_resolution_clock::now(); // Start timing

    if (!m_mesh) {
        std::cerr << "Error: Mesh is not loaded or available." << std::endl;
        return;
    }

    std::cout << "Calculating Symmetry..." << std::endl;

    try {
        MobiusSymmetryDetector detector(m_mesh, m_dijsktra);
        detector.detectSymmetry();
    }
    catch (const std::exception& e) {
        std::cerr << "An error occurred during symmetry detection: " << e.what() << std::endl;
    }

    auto end = std::chrono::high_resolution_clock::now(); // End timing
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "\nSymmetry Calculation Finished." << std::endl;
    std::cout << "Total time for Mobius Symmetry calculation: " << duration.count() << " milliseconds." << std::endl;

    // Append the result to a text file
    std::ofstream outfile("symmetry_timings.txt", std::ios_base::app); // Open in append mode
    if (outfile.is_open()) {
        outfile << std::fixed << std::setprecision(3); // Format to three decimal places
        outfile << "Mobius (" << m_meshName << " "<< m_mesh->samples.size() << " samples) : " << duration.count() / 1000.0 << " seconds" << std::endl;
        outfile.close();
    }
    else {
        std::cerr << "Error: Unable to open timings file for writing." << std::endl;
    }
}
void Display::calculateIntrinsicSymmetry() {
    auto start = std::chrono::high_resolution_clock::now(); // Start timing

    if (!m_mesh) { // Or however you check if m_mesh is valid
        std::cerr << "Error: Mesh is not loaded or available." << std::endl;
        return;
    }

    std::cout << "Calculating Symmetry..." << std::endl;

    // Construct the symmetry detector using the mesh.
    // Assuming 'symm' is the namespace for IntrinsicSymmetryDetector
    symm::IntrinsicSymmetryDetector detector(m_mesh); // If m_mesh is a smart ptr, use .get()
    // If m_mesh is already a raw ptr, just use m_mesh

    // Assemble the Laplacian and mass.
    detector.buildLaplacianAndMass();

    // Compute the eigen-decomposition using (for example) 15 non-trivial eigenpairs.
    // This '15' corresponds to 'd_sig' (signature dimension).
    int signature_dimension = 15;
    if (!detector.computeEigenDecomposition(signature_dimension)) {
        std::cerr << "Failed to compute eigen-decomposition." << std::endl;
        return;
    }

    // Generate candidate sign sequences using the heuristic method.
    // Let's try to get the top 3 promising sign sequences.
    // The first one is usually the identity (all +1s).
    int num_top_sequences_to_find = 2;
    detector.generateSignSequences(num_top_sequences_to_find);
    const auto& sign_sequences = detector.getSignSequences();

    if (sign_sequences.empty()) {
        std::cout << "No sign sequences were generated or selected by the heuristic." << std::endl;
        std::cout << "Symmetry Calculation Finished (no sequences found)." << std::endl;
        return;
    }

    std::cout << "\nFound " << sign_sequences.size() << " candidate sign sequence(s)." << std::endl;

    // --- Process Correspondences for Selected Sign Sequence(s) ---
    // You can choose which sequence(s) to process.
    // For example, to process the first non-identity sequence if available,
    // or iterate through all of them.

    // Let's process correspondences for the first non-identity sequence found by the heuristic.
    // The heuristic implementation adds identity (all +1s) as the first sequence if others are found.
    // So, sign_sequences[0] is likely identity. sign_sequences[1] would be the "best" non-trivial one.

 // --- BEGIN REPLACEMENT in Display::calculateSymmetry ---
// Replace the section that starts with "Let's process correspondences..."

    std::cout << "\nFound " << sign_sequences.size() << " candidate sign sequence(s)." << std::endl;
    std::cout << "Searching for the best non-trivial symmetry..." << std::endl;

    if (sign_sequences.empty()) {
        std::cout << "No sign sequences were found to process." << std::endl;
        return;
    }

    // --- Robustly find the best non-trivial sign sequence ---
    const Eigen::VectorXi* best_S_ptr = nullptr;
    double min_error = std::numeric_limits<double>::infinity();

    // Define the identity sequence for comparison. Assumes all sequences have the same size.
    Eigen::VectorXi identity_S = Eigen::VectorXi::Ones(sign_sequences[0].size());

    for (const auto& current_S : sign_sequences) {
        // Check if the current sequence is the trivial identity sequence
        if (current_S.isApprox(identity_S)) {
            std::cout << "  - Skipping identity sequence." << std::endl;
            continue;
        }

        // For a non-trivial sequence, compute its error score E(S)
        double current_error = detector.computeErrorForSignSequence(current_S);
        std::cout << "  - Evaluating sequence with error: " << current_error << std::endl;

        if (current_error < min_error) {
            min_error = current_error;
            best_S_ptr = &current_S;
        }
    }


    // --- Check if a valid symmetry was found and proceed ---
    if (!best_S_ptr) {
        std::cout << "\nNo non-trivial sign sequences were found." << std::endl;
        std::cout << "Symmetry Calculation Finished." << std::endl;
        return;
    }

    // We now have the best non-trivial sign sequence.
    const Eigen::VectorXi& S_to_use = *best_S_ptr;

    std::cout << "\nProcessing correspondences for best non-trivial Sign Sequence (Error: " << min_error << "):" << std::endl;
    std::cout << "Using Sign Sequence: [";
    for (int i = 0; i < S_to_use.size(); ++i) {
        std::cout << S_to_use(i) << (i == S_to_use.size() - 1 ? "" : " ");
    }
    std::cout << "]" << std::endl;

    std::vector<int> correspondences = detector.computeSymmetryCorrespondences(S_to_use);


    if (correspondences.empty()) {
        std::cout << "Failed to compute correspondences for the selected sign sequence." << std::endl;
    }
    else {
        std::cout << "\nSymmetry correspondences computed (total: " << correspondences.size() << "):" << std::endl;

        // Your logic to store/display some correspondences:
        m_mesh->symmetryAxisVertices.clear(); // Clear previous results
        // --- NEW: Find and store the longest symmetry axis ---
        std::cout << "Finding longest continuous symmetry axis..." << std::endl;
        m_mesh->symmetryAxisVertices = detector.findLongestSymmetryAxis(correspondences);

        if (!m_mesh->symm1.empty()) {
            std::cout << "Found a symmetry axis with " << m_mesh->symmetryAxisVertices.size() << " vertices. These will be colored." << std::endl;
        }
        else {
            std::cout << "No continuous symmetry axis found." << std::endl;
        }
        // --- END NEW SECTION ---
        m_mesh->symm2.clear();

        int step = correspondences.size() / 200;
        int currStep = step;


        std::cout << "Displaying some correspondences (vertex_original -> vertex_symmetric):" << std::endl;
        for (size_t i = 0; i < correspondences.size(); ++i) {
            if (correspondences[i] != -1) { // Check for valid correspondence
                // Your existing logic for specific range
                if (i >= currStep)
                {
                    currStep += step;
                    std::cout << "Vertex " << i << " corresponds to Vertex " << correspondences[i] << std::endl;

                    m_mesh->symm1.push_back(correspondences[i]);
                    m_mesh->symm2.push_back(i);
                }

            }
        }

    }

    auto end = std::chrono::high_resolution_clock::now(); // End timing
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "\nSymmetry Calculation Finished." << std::endl;
    std::cout << "Total time for Intrinsic Symmetry calculation: " << duration.count() << " milliseconds." << std::endl;

    // Append the result to a text file
    std::ofstream outfile("symmetry_timings.txt", std::ios_base::app); // Open in append mode
    if (outfile.is_open()) {
        outfile << std::fixed << std::setprecision(3); // Format to three decimal places
        outfile << "Intrinsic (" << m_meshName << ") : " << duration.count() / 1000.0 << " seconds" << std::endl;
        outfile.close();
    }
    else {
        std::cerr << "Error: Unable to open timings file for writing." << std::endl;
    }
}
