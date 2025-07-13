
#include "Sampling.h"
#include <random>

Sampling::Sampling(Mesh* a_mesh, Dijkstra* a_dijkstra, int a_numOfSamplePts)
	: m_Mesh(a_mesh), m_dijkstra(a_dijkstra), numOfSamplePts(a_numOfSamplePts)
{
	numOfVertices = m_Mesh->vertices.size();
}
void Sampling::calculateSamples()
{  
	std::vector<int> samplePoints = farthestPointSampling();
	m_Mesh->samples = samplePoints;
}

int Sampling::calculateCentroid() 
{
	std::random_device rd;
	std::mt19937 generator(rd());
	std::uniform_int_distribution<int> distribution(0, numOfVertices);

	int initialPoint = distribution(generator);

	std::vector<float> accPathLengths(numOfVertices, 0.0f);
	std::vector<float> minPathLengths(numOfVertices, std::numeric_limits<float>::max());

	std::vector<float> shortestDistance;

	for (int i = 0; i < 100; ++i) {

		shortestDistance = m_dijkstra->getDistanceVector(initialPoint);
		for (int j = 0; j < numOfVertices; ++j) {                  
			accPathLengths[j] += shortestDistance[j];   
			if (shortestDistance[j] < minPathLengths[j])
				minPathLengths[j] = shortestDistance[j];
		}
		initialPoint = globalMaxima(minPathLengths);                   
	}

	return globalMaxima(accPathLengths);
}

std::vector<int> Sampling::farthestPointSampling()
{
	std::vector<int> sampleVertices(numOfSamplePts);
	sampleVertices[0] = calculateCentroid();

	std::vector<float> accMinDist(numOfVertices);
	accMinDist = m_dijkstra->getDistanceVector(sampleVertices[0]);

	sampleVertices[1] = globalMaxima(accMinDist);

	for (int i = 1; i < numOfSamplePts - 1; i++) {

		std::vector<float> shortestDistance = m_dijkstra->getDistanceVector(sampleVertices[i]);

		for (int j = 0; j < numOfVertices; ++j)              
			if (shortestDistance[j] < accMinDist[j])
				accMinDist[j] = shortestDistance[j];

		sampleVertices[i + 1] = globalMaxima(accMinDist); 
	}

	return sampleVertices;
}

int Sampling::globalMaxima(std::vector<float> vals)
{

	if (vals.empty()) 
		return -1;

	float maxVal = vals[0];
	int maxInd = 0;  

	for (int i = 1; i < vals.size(); ++i) {
		if (maxVal < vals[i]) {
			maxInd = i;
			maxVal = vals[i];
		}
	}
	return maxInd;
}
