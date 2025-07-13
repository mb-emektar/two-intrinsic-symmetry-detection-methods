#pragma once
#include "Mesh.h"
#include "Dijkstra.h"

using namespace mesh;

class Sampling
{
public:
	Sampling(Mesh* a_mesh, Dijkstra* a_dijkstra, int a_numOfSamplePts);
	void calculateSamples();
	std::vector<int> farthestPointSampling();
	int calculateCentroid();
	int globalMaxima(std::vector<float> vals);
private:
	Dijkstra* m_dijkstra;
	Mesh* m_Mesh;
	int numOfVertices;
	int numOfSamplePts;
};

