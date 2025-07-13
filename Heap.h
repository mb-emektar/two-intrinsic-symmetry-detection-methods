// C++ Program for implementing Min Heap
#include <iostream>
#include <vector>
#include <climits>

//https://www.geeksforgeeks.org/min-heap-in-cpp/
// Template for representing MinHeap
template <typename T>
class MinHeap {
private:
    std::vector<T> array;
    int size;
    int capacity;

public:
    // Constructor to set the heap's initial capacity
    MinHeap(int capacity) {
        this->size = 0;
        this->capacity = capacity;
        this->array.resize(capacity);
    }

    // Function to restore heap order at index i
    void heapify(int i) {
        // Initialize smallest as root
        int smallest = i;
        //  Find the Left child index
        int left = 2 * i + 1;
        // Find the Right child index
        int right = 2 * i + 2;

        // If left child is smaller than root
        if (left < size && array[left] < array[smallest])
            smallest = left;

        // If right child is smaller than the smallest so far
        if (right < size && array[right] < array[smallest])
            smallest = right;

        // If smallest is not root
        if (smallest != i) {
            std::swap(array[i], array[smallest]);
            heapify(smallest);
        }
    }

    // Function to create a min heap from a given array
    void buildHeap(const std::vector<T>& arr) {
        capacity = arr.size();
        size = capacity;
        array = arr;

        // Heapify the (n-1)/2 nodes
        for (int i = (size - 1) / 2; i >= 0; i--) {
            heapify(i);
        }
    }

    // Function to insert a new node into the min heap
    void insertNode(T value) {
        if (size == capacity) {
            // Resize the heap if necessary
            capacity *= 2;
            array.resize(capacity);
        }

        size++;
        int i = size - 1;
        array[i] = value;

        while (i != 0 && array[(i - 1) / 2] > array[i]) {
            std::swap(array[i], array[(i - 1) / 2]);
            i = (i - 1) / 2;
        }
    }

    // Function to get the topmost value from the min heap
    T peek() {
        if (size <= 0)
            return -1;  // Indicates that the heap is empty

        return array[0];
    }

    // Function to remove and return the minimum value from the heap
    T extractMin() {
        if (size <= 0)
            return -1;  // Indicates that the heap is empty
        if (size == 1) {
            size--;
            return array[0];
        }

        // Store the minimum value, and remove it
        T root = array[0];
        array[0] = array[size - 1];
        size--;
        // Heapify the root node after deletion
        heapify(0);
        return root;
    }

    // Function to delete a specific node from the heap
    void DeleteNode(T key) {
        // Find the index of the key
        int index = -1;
        for (int i = 0; i < size; ++i) {
            if (array[i] == key) {
                index = i;
                break;
            }
        }
        // If key is not found, return
        if (index == -1) {
            std::cout << "Key not found" << std::endl;
            return;
        }

        if (index == size - 1) {
            size--;
            return;
        }

        // Move the last element to the index of the key to be deleted
        array[index] = array[size - 1];
        size--;

        // Heapify down to heapify the rest of the min heap
        heapify(index);
    }

    // Function to print the values of  the  min heap
    void printHeap() const {
        for (int i = 0; i < size; ++i)
            std::cout << array[i] << " ";
        std::cout << std::endl;
    }
};

