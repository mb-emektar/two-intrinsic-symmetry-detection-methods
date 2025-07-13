#include "Mesh.h"

namespace mesh {

Mesh::Mesh()
{

}
Mesh::Mesh(char* name) : m_name(name)
{
	FILE* fPtr = fopen(name, "r");
	char str[334];

	fscanf(fPtr, "%s", str);

	int nVerts, nTris, n, i = 0;
	float x, y, z;

	fscanf(fPtr, "%d %d %d\n", &nVerts, &nTris, &n);
	while (i++ < nVerts)
	{
		fscanf(fPtr, "%f %f %f", &x, &y, &z);
		addVertex(x, y, z);
	}

	while (fscanf(fPtr, "%d", &i) != EOF)
	{
		fscanf(fPtr, "%f %f %f", &x, &y, &z);
		addTriangle((int)x, (int)y, (int)z);
	}

	fclose(fPtr);
}

Mesh::~Mesh()
{
}
Mesh::Mesh(Mesh& other) {
	vertices = other.vertices;
	tris = other.tris;
	edges = other.edges;
	highlightedEdges = other.highlightedEdges;
	samples = other.samples;
}
void Mesh::saveAsFile() {
	FILE* fPtr = fopen("outputs/Output.off", "w");
	fprintf(fPtr, "OFF\n");
	fprintf(fPtr, "%d %d %d\n", vertices.size(), tris.size(), edges.size());
	for (int i = 0; i < vertices.size(); i++)
		fprintf(fPtr, "%f %f %f\n", vertices[i]->coords.x, vertices[i]->coords.y, vertices[i]->coords.z);
	for (int i = 0; i < tris.size(); i++)
		fprintf(fPtr, "3 %d %d %d\n", tris[i]->v1i, tris[i]->v2i, tris[i]->v3i);
	fclose(fPtr);
}
float Mesh::getTotalArea() {
	float total_area = 0.0f;
	for (const auto& tri : tris) {
		const Vec3& v1 = vertices[tri->v1i]->coords;
		const Vec3& v2 = vertices[tri->v2i]->coords;
		const Vec3& v3 = vertices[tri->v3i]->coords;
		total_area += 0.5f * (v2 - v1).cross(v3 - v1).length();
	}
	return total_area;

}
void Mesh::createCube(float sideLen)
{
	//coordinates
	float flbc[3] = { 0, 0, 0 }, deltaX = 0, deltaY = 0, deltaZ = 0;
	for (int v = 0; v < 8; v++)
	{
		switch (v)
		{
		case 1:
			deltaX = sideLen;
			break;
		case 2:
			deltaZ = -sideLen;
			break;
		case 3:
			deltaX = 0;
			break;
		case 4:
			deltaZ = 0;
			deltaY = sideLen;
			break;
		case 5:
			deltaX = sideLen;
			break;
		case 6:
			deltaZ = -sideLen;
			break;
		default:
			deltaX = 0;;
			break;
		}
		addVertex(flbc[0] + deltaX, flbc[1] + deltaY, flbc[2] + deltaZ);
	}

	addTriangle(0, 2, 1);
	addTriangle(0, 3, 2);

	addTriangle(1, 2, 5);
	addTriangle(2, 6, 5);

	addTriangle(2, 3, 6);
	addTriangle(3, 7, 6);

	addTriangle(3, 4, 7);
	addTriangle(3, 0, 4);

	addTriangle(4, 5, 6);
	addTriangle(4, 6, 7);

	addTriangle(0, 1, 5);
	addTriangle(0, 5, 4);
}

void Mesh::addTriangle(int v1, int v2, int v3)
{
	int idx = tris.size();
	tris.push_back(new Triangle(idx, v1, v2, v3));

	//set up structure

	vertices[v1]->triList.push_back(idx);
	vertices[v2]->triList.push_back(idx);
	vertices[v3]->triList.push_back(idx);

	if (!makeVertsNeighbor(v1, v2))
		addEdge(v1, v2);

	if (!makeVertsNeighbor(v1, v3))
		addEdge(v1, v3);

	if (!makeVertsNeighbor(v2, v3))
		addEdge(v2, v3);

}

bool Mesh::makeVertsNeighbor(int v1i, int v2i)
{
	//returns true if v1i already neighbor w/ v2i; false o/w

	for (int i = 0; i < vertices[v1i]->vertList.size(); i++)
		if (vertices[v1i]->vertList[i] == v2i)
			return true;


	vertices[v1i]->vertList.push_back(v2i);
	vertices[v2i]->vertList.push_back(v1i);
	return false;
}

void Mesh::addVertex(float x, float y, float z)
{
	int idx = vertices.size();
	//	float* c = new float[3];
	//	c[0] = x;
	//	c[1] = y;
	//	c[2] = z;

	vertices.push_back(new Vertex(idx, x, y, z));
}

void Mesh::addEdge(int v1, int v2)
{
	int idx = edges.size();

	edges.push_back(new Edge(idx, v1, v2));

	vertices[v1]->edgeList.push_back(idx);
	vertices[v2]->edgeList.push_back(idx);
}

float Mesh::getLength(int v1, int v2) {
	return
		sqrt(pow((vertices[v2]->coords.x - vertices[v1]->coords.x), 2) +
			pow((vertices[v2]->coords.y - vertices[v1]->coords.y), 2) +
			pow((vertices[v2]->coords.z - vertices[v1]->coords.z), 2));
}

void Mesh::normalizeMesh() {
	// 1. Check if there are any vertices to process
	if (vertices.empty()) {
		std::cout << "Cannot normalize an empty mesh." << std::endl;
		return;
	}

	// 2. Find the bounding box of the mesh
	Vec3 min_coords = vertices[0]->coords;
	Vec3 max_coords = vertices[0]->coords;

	for (size_t i = 1; i < vertices.size(); ++i) {
		const Vec3& p = vertices[i]->coords;
		min_coords.x = std::min(min_coords.x, p.x);
		min_coords.y = std::min(min_coords.y, p.y);
		min_coords.z = std::min(min_coords.z, p.z);
		max_coords.x = std::max(max_coords.x, p.x);
		max_coords.y = std::max(max_coords.y, p.y);
		max_coords.z = std::max(max_coords.z, p.z);
	}

	// 3. Compute the geometric center of the bounding box
	Vec3 center = (min_coords + max_coords) * 0.5f;

	// 4. Translate all vertices to center the mesh at the origin
	//    and find the maximum squared distance from the new origin
	float max_dist_sq = 0.0f;
	for (size_t i = 0; i < vertices.size(); ++i) {
		vertices[i]->coords -= center;

		// Use dot product for squared length, which is more efficient than length()
		float dist_sq = vertices[i]->coords.dot(vertices[i]->coords);
		if (dist_sq > max_dist_sq) {
			max_dist_sq = dist_sq;
		}
	}

	// 5. If max_dist_sq is zero, all points are the same. Avoid division by zero.
	if (max_dist_sq == 0.0f) {
		std::cout << "Mesh has zero extent; cannot scale." << std::endl;
		return;
	}

	// 6. Compute the scaling factor (1 / max_distance)
	float scale_factor = 100.0f / std::sqrt(max_dist_sq);

	// 7. Scale all vertices to fit within a unit sphere
	for (size_t i = 0; i < vertices.size(); ++i) {
		vertices[i]->coords *= scale_factor;
	}

	std::cout << "Mesh has been normalized to fit within a unit sphere." << std::endl;
}

} // namespace mesh