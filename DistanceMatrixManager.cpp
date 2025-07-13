#include "DistanceMatrixManager.h"

DistanceMatrixManager::DistanceMatrixManager(const std::string& meshName, int numVertices)
    : m_num_vertices(numVertices) {
    // Create a unique directory for cache
    m_cache_directory = std::filesystem::path("cache") / meshName;
}

std::string DistanceMatrixManager::getFilePath(int sourceIndex) const {
    return (m_cache_directory / (std::to_string(sourceIndex) + ".dist")).string();
}

bool DistanceMatrixManager::cacheDirectoryExists() const {
    return std::filesystem::exists(m_cache_directory);
}

void DistanceMatrixManager::createCacheDirectory() {
    std::filesystem::create_directories(m_cache_directory);
}

void DistanceMatrixManager::saveVector(int sourceIndex, const std::vector<float>& distances) {
    if (distances.size() != m_num_vertices) {
        std::cerr << "Error: Mismatch in distance vector size during save for index " << sourceIndex << std::endl;
        return;
    }

    std::ofstream outfile(getFilePath(sourceIndex), std::ios::binary);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file for writing: " << getFilePath(sourceIndex) << std::endl;
        return;
    }

    // --- NEW: Write a header ---
    outfile.write(reinterpret_cast<const char*>(&m_num_vertices), sizeof(int));

    outfile.write(reinterpret_cast<const char*>(distances.data()), m_num_vertices * sizeof(float));
    outfile.close();
}

bool DistanceMatrixManager::loadVector(int sourceIndex, std::vector<float>& outDistances) {
    std::ifstream infile(getFilePath(sourceIndex), std::ios::binary);
    if (!infile.is_open()) {
        return false;
    }

    // --- NEW: Read and validate the header ---
    int num_vertices_in_file = 0;
    infile.read(reinterpret_cast<char*>(&num_vertices_in_file), sizeof(int));

    if (!infile || num_vertices_in_file != m_num_vertices) {
        std::cerr << "Warning: Cache file validation failed for index " << sourceIndex
            << ". Expected " << m_num_vertices << " vertices, but file reports "
            << num_vertices_in_file << "." << std::endl;
        infile.close();
        return false;
    }

    outDistances.resize(m_num_vertices);
    infile.read(reinterpret_cast<char*>(outDistances.data()), m_num_vertices * sizeof(float));

    bool read_success = infile.gcount() == m_num_vertices * sizeof(float);
    infile.close();
    return read_success;
}

bool DistanceMatrixManager::vectorExists(int sourceIndex) const {
    return std::filesystem::exists(getFilePath(sourceIndex));
}