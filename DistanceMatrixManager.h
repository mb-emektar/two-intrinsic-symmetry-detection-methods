#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <filesystem>
#include <iostream>

class DistanceMatrixManager {
public:
    DistanceMatrixManager(const std::string& meshName, int numVertices);

    // Saves a single distance vector for a requested vertex.
    void saveVector(int sourceIndex, const std::vector<float>& distances);

    // Loads a single distance vector for a requested vertex.
    bool loadVector(int sourceIndex, std::vector<float>& outDistances);

    // Checks if the pre-computed data for a given vertex exists.
    bool vectorExists(int sourceIndex) const;

    // Checks if the base directory for caching exists.
    bool cacheDirectoryExists() const;

    // Creates the directory for caching.
    void createCacheDirectory();

private:
    std::string getFilePath(int sourceIndex) const;
    std::filesystem::path m_cache_directory;
    int m_num_vertices;
};