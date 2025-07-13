#pragma once

// Standard library includes
#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <list>        
#include <unordered_map>
#include <mutex>     
#include <chrono>
#include <fstream>

// Project-specific includes
#include "Mesh.h"
#include "fiboheap.h"  
#include "DistanceMatrixManager.h" 

using namespace mesh;

/**
 * @class Node
 * @brief A helper struct  used in algorithm for the priority queues
 *
 * This class holds the index of a vertex and its distance from a source.
 * It includes comparison operators, that are necessary for the priority queue
 * to order its elements correctly.
 */
class Node {
public:
    int vertex;
    float distance;

    Node(int a_vertex = -1, float a_distance = -1) : vertex(a_vertex), distance(a_distance) {}

    // Comparison operators for sorting within a priority queue (heap).
    bool operator<(const Node& rhs) const { return distance < rhs.distance; }
    bool operator>(const Node& rhs) const { return distance > rhs.distance; }
    bool operator==(const Node& rhs) const { return vertex == rhs.vertex; }
    bool operator==(int rhs) const { return vertex == rhs; }
};


/**
 * @class Dijkstra
 * @brief Implements Dijkstra's algorithm to determine shortest paths for all vertices.
 *
 * This class offers fibo-heap implementation of Dijkstra's algorithm and features
 * a sophisticated two-level caching method to handle expensive distance calculations efficiently.
 *
 * - Pre-computation step can calculate all-pairs shortest paths and
 * store every distance vector in its own file on filesystem.
 * - Thread-safe in-memory LRU (Least Recently Used) cache provides fast
 * access to frequently used distance vectors in runtime.
 */
class Dijkstra {
public:
    /**
     * @brief Constructor
     * @param a_mesh Pointer to the mesh on which to run the dijkstra.
     */
    Dijkstra(Mesh* a_mesh);

    // --- Core Public Interface ---

    /**
     * @brief Performs an expensive pre-computation of all-pairs shortest paths.
     * This function calculates the distance vector for all vertices and saves each
     * to a file on disk. It is multi-threaded to be efficient. If the cache
     * already exists in directory, this step is skipped.
     */
    void precomputeAndCacheAllDistances();

    /**
     * @brief Brings the distance vector of given vertex.
     * This is the main, thread-safe method for accessing distance vectors. It uses
     * a two-level cache (RAM -> Disk) to provide efficiency.
     * If a vector doesn't exist in any cache, it is computed and saved.
     * @param index The index of the vertex.
     * @return Distance matrix
     */
    std::vector<float> getDistanceVector(int index);

    /**
     * @brief Get access to the mesh.
     * @return Pointer to the Mesh object.
     */
    Mesh* getMesh();

    // --- Different Algorithm Implementations ---

    /**
     * @brief Calculates distance vector using a Fibo heap.
     * @param startingIdx The index of the source vertex.
     * @return Distance vector of the given vertex
     */
    std::vector<float> calculateDijkstraFiboHeap(int startingIdx);

private:
    Mesh* m_Mesh;
    int startNode, endNode;
    int numOfVertices;

    // Manages distance vectors efficienty.
    DistanceMatrixManager* m_matrix_manager;

    // --- Thread-Safe LRU Cache Members ---
    std::mutex m_cache_mutex; // Protects access to RAM cache.

    // Number of distance vectors to keep in RAM.
    const size_t MAX_RAM_CACHE_SIZE = 4000;

    // In-memory cache. Maps a vertex to its distance vector and an
    // iterator pointing to its location in the usage list.
    std::unordered_map<int, std::pair<std::vector<float>, std::list<int>::iterator>> m_distance_ram_cache;

    // Tracks the usage order of keys. Most recently used distance vector is at the front.
    std::list<int> m_cache_usage_order;
};