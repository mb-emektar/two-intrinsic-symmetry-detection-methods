#include "Dijkstra.h"
#include <map>
#include <atomic>
#include <filesystem>
#include <sstream>

// Constructor
Dijkstra::Dijkstra(Mesh* a_mesh) : m_Mesh(a_mesh) {
    numOfVertices = m_Mesh->vertices.size();

    // Extract the base name of mesh file to use for the cache directory.
    std::string mesh_name = std::filesystem::path(m_Mesh->m_name).stem().string();
    m_matrix_manager = new DistanceMatrixManager(mesh_name, numOfVertices);

    startNode = -1;
    endNode = -1;
}

// Use a Fibonacci heap for the priority queue for the most efficient implementation, .
std::vector<float> Dijkstra::calculateDijkstraFiboHeap(int startingIdx) {
    FibonacciHeap<Node> toExplore;
    node<Node>** nodePtr = new node<Node>*[numOfVertices];
    std::vector<float> dist(numOfVertices, std::numeric_limits<float>::max());

    dist[startingIdx] = 0;
    for (int i = 0; i < numOfVertices; ++i) {
        nodePtr[i] = toExplore.insert(Node(i, dist[i]));
    }

    while (!toExplore.isEmpty()) {
        Node u = toExplore.removeMinimum();

        for (auto neighbor : m_Mesh->vertices[u.vertex]->vertList) {
            // Assuming uniform edge weights of 1.
            float newDist = dist[u.vertex]  +m_Mesh->getLength(u.vertex, neighbor);//+ 1;
            if (newDist < dist[neighbor]) {
                dist[neighbor] = newDist;
                toExplore.decreaseKey(nodePtr[neighbor], Node(neighbor, newDist));
            }
        }
    }

    delete[] nodePtr;
    return dist;
}

// Retrieves a distance vector with caching method.
std::vector<float> Dijkstra::getDistanceVector(int index) {
    // lock_guard ensures the mutex is locked automatically
    // unlocked when the function returns
    std::lock_guard<std::mutex> lock(m_cache_mutex);

    // 1. Check RAM cache first (LRU cache).
    auto it = m_distance_ram_cache.find(index);
    if (it != m_distance_ram_cache.end()) {
        m_cache_usage_order.splice(m_cache_usage_order.begin(), m_cache_usage_order, it->second.second);
        return it->second.first;
    }

    // --- RAM Cache Miss ---

    // 2. Try to load the vector from the disk cache.
    std::vector<float> distanceVector;
    if (m_matrix_manager->loadVector(index, distanceVector)) {
        if (m_distance_ram_cache.size() >= MAX_RAM_CACHE_SIZE) {
            int key_to_evict = m_cache_usage_order.back();
            m_cache_usage_order.pop_back();
            m_distance_ram_cache.erase(key_to_evict);
        }
        m_cache_usage_order.push_front(index);
        m_distance_ram_cache[index] = { distanceVector, m_cache_usage_order.begin() };
        return distanceVector;
    }
    
    // --- DISK Cache Miss ---

    // 3. Complete Cache Miss: The vector is not in RAM or on disk.
    // Calculate and save it to disk and RAM for future use.
    std::cout << "Warning: Vector " << index << " not found in cache. Calculating on-the-fly." << std::endl;
    distanceVector = calculateDijkstraFiboHeap(index);
    m_matrix_manager->saveVector(index, distanceVector); // save to disk.

    // add to RAM cache
    if (m_distance_ram_cache.size() >= MAX_RAM_CACHE_SIZE) {
        int key_to_evict = m_cache_usage_order.back();
        m_cache_usage_order.pop_back();
        m_distance_ram_cache.erase(key_to_evict);
    }
    m_cache_usage_order.push_front(index);
    m_distance_ram_cache[index] = { distanceVector, m_cache_usage_order.begin() };

    return distanceVector;
}

// Returns a pointer to the mesh object.
Mesh* Dijkstra::getMesh() {
    return m_Mesh;
}

// Performs the one-time, expensive pre-computation of all distance vectors,
// saving them to the disk cache for future runs.
void Dijkstra::precomputeAndCacheAllDistances() {
    if (m_matrix_manager->cacheDirectoryExists()) {
        std::cout << "Distance matrix cache found. Skipping pre-computation." << std::endl;
        return;
    }

    std::cout << "No distance cache found. Starting pre-computation..." << std::endl;
    m_matrix_manager->createCacheDirectory();

    std::atomic<int> progress_counter(0);

    // Use OpenMP to parallelize the loop
#pragma omp parallel for
    for (int i = 0; i < numOfVertices; ++i) {
        std::vector<float> distanceVector = calculateDijkstraFiboHeap(i);
        m_matrix_manager->saveVector(i, distanceVector); // Save each vector to its own file.

        int completed = ++progress_counter;
        if (completed % 100 == 0 || completed == numOfVertices) {
#pragma omp critical
            {
                std::cout << "Pre-computing distances: " << (completed * 100) / numOfVertices << "%\r" << std::flush;
            }
        }
    }
    std::cout << std::endl << "Pre-computation complete. Matrix cached to disk." << std::endl;
}
