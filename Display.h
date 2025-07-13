

#pragma once
#include "IntrinsicSymmetryDetector.h"
#include "MobiusSymmetryDetector.h"
#define _CRT_SECURE_NO_DEPRECATE
#include <inttypes.h>
#define HAVE_INT8_T
#include <cmath>
#include <xlocnum>


#include <Inventor/Win/SoWin.h>
#include <Inventor/Win/viewers/SoWinExaminerViewer.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoCone.h>

#include "Mesh.h"
#include "Dijkstra.h"
#include "Painter.h"
#include "Sampling.h"
#include "DistanceMatrixManager.h"


using namespace symm;
using namespace mesh;

class Display
{
public:

	Display(char* argv,char* name);
	~Display();
	void displayMesh();
	void createMembers(int methodSelection);
	void calculateSampling();

	void calculateMobiusSymmetry();
	void calculateIntrinsicSymmetry();

private:
	char* m_meshName;
	HWND window;
	SoWinExaminerViewer* viewer;
	SoSeparator* root;

	Sampling* m_sampling;
	Dijkstra* m_dijsktra;
	Mesh* m_mesh;
	Painter* m_painter;

	int numOfSamplePts;
	std::vector<std::vector<float>> m_distanceMatrix;
};

