#include "Painter.h"
#include "Dijkstra.h"

Painter::Painter() 
{

	for (int i = 1; i < 8; ++i) {

		int a = (i >> 2) & 1; // Get the bit at the 2nd position (from the right, 0-indexed)
		int c = i & 1;       // Get the bit at the 0th position
		int b = 0;
		if(a || c)
			b = (i >> 1) & 1; // Get the bit at the 1st position
		SbVec3f color(a, b, c);
		colors.push_back(color);
	}
}

SoSeparator* Painter::getSphereSep(Mesh* mesh, float deltaX, float deltaY, float scale)
{
	//returns a set of spheres to highlight each mesh.samples[i]

	SoSeparator* spheresSep = new SoSeparator();

	float radius = RADIUS;

	for (int i = 0; i < (int) mesh->samples.size(); i++)
	{
		//1 sphere for this sample
		SoSeparator* sphere1Sep = new SoSeparator;

		//transformation
		SoTransform* tra = new SoTransform();
		float * temp = mesh->vertices[ mesh->samples[i] ]->coords.toFloat3();
		tra->translation.setValue(scale*temp[0]+deltaX, scale*temp[1]+deltaY, scale*temp[2]);
		sphere1Sep->addChild(tra);
		delete[] temp;

		//material
		SoMaterial* ma = new SoMaterial;
		if (i == 0)
			ma->diffuseColor.setValue(SbColor(0.0f, 0.0f, 0.7f));
		else if (i == 1)
			ma->diffuseColor.setValue(SbColor(0.0f, 0.0f, 0.0f));
		else if (i == 2)
			ma->diffuseColor.setValue(SbColor(0.0f, 0.7f, 0.0f));
		else if (i == 3)
			ma->diffuseColor.setValue(SbColor(0.7f, 0.0f, 0.7f));
		else if (i == 4)
			ma->diffuseColor.setValue(SbColor(0.7f, 0.7f, 0.0f));
		else
			ma->diffuseColor.setValue(SbColor(0.7f, 0.0f, 0.0f));

		sphere1Sep->addChild(ma);
		
		//shape
		SoSphere* sph1 = new SoSphere();
		sph1->radius = radius;
		sphere1Sep->addChild(sph1); //whose position is decided by the translation applied above

		spheresSep->addChild(sphere1Sep);
	}
	
	return spheresSep;
}



SoSeparator* Painter::getShapeSep(Mesh* mesh)
{
	SoSeparator* res = new SoSeparator();

	//transformation
	//not needed

	//color
	SoMaterial* mat = new SoMaterial();
	mat->diffuseColor.setValue(0, 1, 1); //paint all vertices with this color
	//mat->transparency = 0.5f : 0.0f; //0 makes it completely opaque, the default

	bool youWantToPaintEachVertexDifferently = false;
	if (youWantToPaintEachVertexDifferently)
		for (int i = 0; i < (int) mesh->vertices.size(); i++) //i = 0 obj->color above overwritten here
			mat->diffuseColor.set1Value(i, mesh->vertices[i]->color); //vert color according to its x-y-z coord (for mesh1) and to the transferred color (for mesh2)


	res->addChild(mat);

	SoShapeHints* hints = new SoShapeHints;
	hints->creaseAngle = 3.14;
	res->addChild(hints); //Gouraud shading

	if (youWantToPaintEachVertexDifferently)
	{
		SoMaterialBinding* materialBinding = new SoMaterialBinding; //for 2+ diffuse color usage on the same mesh
		materialBinding->value = SoMaterialBinding::PER_VERTEX_INDEXED;
		res->addChild(materialBinding);
	}

	//shape
	SoCoordinate3* coords = new SoCoordinate3();
	for (int c = 0; c < mesh->vertices.size(); c++)
		coords->point.set1Value(c, mesh->vertices[c]->coords[0], mesh->vertices[c]->coords[1], mesh->vertices[c]->coords[2]);
	SoIndexedFaceSet* faceSet = new SoIndexedFaceSet();
	for (int c = 0; c < mesh->tris.size(); c++)
	{
		faceSet->coordIndex.set1Value(c*4, mesh->tris[c]->v1i);
		faceSet->coordIndex.set1Value(c*4 + 1, mesh->tris[c]->v2i);
		faceSet->coordIndex.set1Value(c*4 + 2, mesh->tris[c]->v3i);
		faceSet->coordIndex.set1Value(c*4 + 3, -1);

		if (youWantToPaintEachVertexDifferently)
		{
			faceSet->materialIndex.set1Value(0 + 4*c, mesh->tris[c]->v1i);
			faceSet->materialIndex.set1Value(1 + 4*c, mesh->tris[c]->v2i);
			faceSet->materialIndex.set1Value(2 + 4*c, mesh->tris[c]->v3i);
		}
	}
	res->addChild(coords);
	res->addChild(faceSet);



	for (int k = 0; k < mesh->highlightedEdges.size(); k++)
	{
		bool drawThickEdges = false;
		if (mesh->highlightedEdges[k].size())
			drawThickEdges = true;

		if (drawThickEdges) //draw thick edges (may be useful in geodesic path drawing)
		{
			float deltaX = 0;

			SoSeparator* thickEdgeSep = new SoSeparator;
			//material
			SoMaterial* ma = new SoMaterial;
			SbVec3f currColor = colors[k];
			ma->diffuseColor.set1Value(0, currColor[2], currColor[1], currColor[0]);
			thickEdgeSep->addChild(ma);
			SoDrawStyle* sty = new SoDrawStyle;	sty->lineWidth = 5.0f;	thickEdgeSep->addChild(sty);

			//shape
			SoIndexedLineSet* ils = new SoIndexedLineSet;
			SoCoordinate3* co = new SoCoordinate3;

			//assumes no edge in sedges is removed

			//		std::cout << "" << std::endl;

			for (unsigned int se = 0; se < mesh->highlightedEdges[k].size(); se++)
			{
				float* temp1 = mesh->vertices[mesh->highlightedEdges[k][se]->v1i]->coords.toFloat3(),
					* temp2 = mesh->vertices[mesh->highlightedEdges[k][se]->v2i]->coords.toFloat3();

				SbVec3f end1 = temp1 + SbVec3f(deltaX, 0.0f, 0.0f),
					end2 = temp2 + SbVec3f(deltaX, 0.0f, 0.0f);

				co->point.set1Value(2 * se, end1);
				co->point.set1Value(2 * se + 1, end2);
				delete[] temp1;
				delete[] temp2;
			}
			//		std::cout << "" << std::endl;

			for (unsigned int ci = 0; ci < mesh->highlightedEdges[k].size(); ci++)
			{
				ils->coordIndex.set1Value(3 * ci, 2 * ci);	ils->coordIndex.set1Value(3 * ci + 1, 2 * ci + 1);
				ils->coordIndex.set1Value(3 * ci + 2, -1); //end this edge with -1
			}
			thickEdgeSep->addChild(co);	thickEdgeSep->addChild(ils);
			res->addChild(thickEdgeSep);
		}
	}

	// Highlight vertices
	SoSeparator* highlightedVertexSep = new SoSeparator();
	{	
		// Material for highlighted vertices
		SoMaterial* highlightMat = new SoMaterial();
		highlightMat->diffuseColor.setValue(0, 1, 1); // Highlight color (e.g., green)
		highlightedVertexSep->addChild(highlightMat);

		// Point size for highlighted vertices
		SoDrawStyle* vertexStyle = new SoDrawStyle();
		vertexStyle->pointSize = 8.0f; // Adjust size as needed
		highlightedVertexSep->addChild(vertexStyle);

		// Coordinates for highlighted vertices
		SoCoordinate3* highlightCoords = new SoCoordinate3();
		for (int hv = 0; hv < mesh->samples.size(); ++hv)
		{
			highlightCoords->point.set1Value(hv,
				mesh->vertices[mesh->samples[hv]]->coords[0],
				mesh->vertices[mesh->samples[hv]]->coords[1],
				mesh->vertices[mesh->samples[hv]]->coords[2]);
		}
		highlightedVertexSep->addChild(highlightCoords);

		// Shape for highlighted vertices
		SoPointSet* pointSet = new SoPointSet();
		highlightedVertexSep->addChild(pointSet);
	}


	res->addChild(highlightedVertexSep);
	///////////////////
	/*
	// Highlight vertices
	highlightedVertexSep = nullptr;
	highlightedVertexSep = new SoSeparator();
	{
		// Material for highlighted vertices
		SoMaterial* highlightMat = new SoMaterial();
		highlightMat->diffuseColor.setValue(1, 1, 1);
		highlightedVertexSep->addChild(highlightMat);

		// Point size for highlighted vertices
		SoDrawStyle* vertexStyle = new SoDrawStyle();
		vertexStyle->pointSize = 8.0f; // Adjust size as needed
		highlightedVertexSep->addChild(vertexStyle);

		// Coordinates for highlighted vertices
		SoCoordinate3* highlightCoords = new SoCoordinate3();
		for (int hv = 0; hv < mesh->symm1.size(); ++hv)
		{
			highlightCoords->point.set1Value(hv,
				mesh->vertices[mesh->symm1[hv]]->coords[0],
				mesh->vertices[mesh->symm1[hv]]->coords[1],
				mesh->vertices[mesh->symm1[hv]]->coords[2]);
		}
		highlightedVertexSep->addChild(highlightCoords);

		// Shape for highlighted vertices
		SoPointSet* pointSet = new SoPointSet();
		highlightedVertexSep->addChild(pointSet);
	}

	res->addChild(highlightedVertexSep);

	// Highlight vertices
	highlightedVertexSep = nullptr;
	highlightedVertexSep = new SoSeparator();
	{
		// Material for highlighted vertices
		SoMaterial* highlightMat = new SoMaterial();
		highlightMat->diffuseColor.setValue(0, 0, 0);
		highlightedVertexSep->addChild(highlightMat);

		// Point size for highlighted vertices
		SoDrawStyle* vertexStyle = new SoDrawStyle();
		vertexStyle->pointSize = 8.0f; // Adjust size as needed
		highlightedVertexSep->addChild(vertexStyle);

		// Coordinates for highlighted vertices
		SoCoordinate3* highlightCoords = new SoCoordinate3();
		for (int hv = 0; hv < mesh->symm2.size(); ++hv)
		{
			highlightCoords->point.set1Value(hv,
				mesh->vertices[mesh->symm2[hv]]->coords[0],
				mesh->vertices[mesh->symm2[hv]]->coords[1],
				mesh->vertices[mesh->symm2[hv]]->coords[2]);
		}
		highlightedVertexSep->addChild(highlightCoords);

		// Shape for highlighted vertices
		SoPointSet* pointSet = new SoPointSet();
		highlightedVertexSep->addChild(pointSet);
	}

	res->addChild(highlightedVertexSep);
	*/
	// Highlight vertices
	highlightedVertexSep = nullptr;
	highlightedVertexSep = new SoSeparator();
	{
		// Material for highlighted vertices
		SoMaterial* highlightMat = new SoMaterial();
		highlightMat->diffuseColor.setValue(1, 0, 0);
		highlightedVertexSep->addChild(highlightMat);

		// Point size for highlighted vertices
		SoDrawStyle* vertexStyle = new SoDrawStyle();
		vertexStyle->pointSize = 8.0f; // Adjust size as needed
		highlightedVertexSep->addChild(vertexStyle);

		// Coordinates for highlighted vertices
		SoCoordinate3* highlightCoords = new SoCoordinate3();
		for (int hv = 0; hv < mesh->symmetryAxisVertices.size(); ++hv)
		{
			highlightCoords->point.set1Value(hv,
				mesh->vertices[mesh->symmetryAxisVertices[hv]]->coords[0],
				mesh->vertices[mesh->symmetryAxisVertices[hv]]->coords[1],
				mesh->vertices[mesh->symmetryAxisVertices[hv]]->coords[2]);
		}
		highlightedVertexSep->addChild(highlightCoords);

		// Shape for highlighted vertices
		SoPointSet* pointSet = new SoPointSet();
		highlightedVertexSep->addChild(pointSet);
	}

	res->addChild(highlightedVertexSep);

	// MODIFICATION STARTS HERE: Draw lines between symm1 and symm2
// Remove the old SoPointSet drawing for symm1 and symm2

	SoSeparator* symmLinesSep = new SoSeparator();
	{
		// Material for the lines
		SoMaterial* lineMat = new SoMaterial();
		lineMat->diffuseColor.setValue(1.0f, 1.0f, 0.0f); // Yellow color for lines
		symmLinesSep->addChild(lineMat);

		// Draw style for line width
		SoDrawStyle* lineStyle = new SoDrawStyle();
		lineStyle->lineWidth = 2.0f; // Set line width
		symmLinesSep->addChild(lineStyle);

		SoCoordinate3* lineCoords = new SoCoordinate3();
		SoIndexedLineSet* lineSet = new SoIndexedLineSet();

		int num_symm1_pts = static_cast<int>(mesh->symm1.size());
		int num_symm2_pts = static_cast<int>(mesh->symm2.size());
		int num_lines_to_draw = std::min(num_symm1_pts, num_symm2_pts);

		if (num_symm1_pts != num_symm2_pts) {
			std::cerr << "Painter::getShapeSep: Warning! symm1 (" << num_symm1_pts
				<< ") and symm2 (" << num_symm2_pts
				<< ") have different sizes. Drawing lines for the minimum count: "
				<< num_lines_to_draw << std::endl;
		}

		if (num_lines_to_draw > 0) {
			for (int i = 0; i < num_lines_to_draw; ++i) {
				// Assuming mesh->vertices[index]->coords gives an SbVec3f or similar
				// and symm1/symm2 contain valid indices for mesh->vertices
				Vec3 p1 = mesh->vertices[mesh->symm1[i]]->coords;
				Vec3 p2 = mesh->vertices[mesh->symm2[i]]->coords;

				lineCoords->point.set1Value(2 * i, p1[0],p1[1],	p1[2]);

				lineCoords->point.set1Value(2 * i + 1, p2[0], p2[1], p2[2]);

				lineSet->coordIndex.set1Value(3 * i, 2 * i);       // Start index of line in lineCoords
				lineSet->coordIndex.set1Value(3 * i + 1, 2 * i + 1); // End index of line in lineCoords
				lineSet->coordIndex.set1Value(3 * i + 2, -1);      // End of this line segment
			}
			symmLinesSep->addChild(lineCoords);
			symmLinesSep->addChild(lineSet);
		}
		// If num_lines_to_draw is 0, an empty separator is added, which is fine.
	}
	res->addChild(symmLinesSep); // Add the separator for symm lines to the result
	// MODIFICATION ENDS HERE
	///////////////////
	// 
	
	/*
	// Text for vertex indices
	SoSeparator* indexTextSep = new SoSeparator();
	SoFont* indexFont = new SoFont();
	indexFont->size.setValue(40.0f); // Adjust font size as needed
	indexTextSep->addChild(indexFont);

	// Material for text color
	SoMaterial* textColorMat = new SoMaterial();
	textColorMat->diffuseColor.setValue(1, 1, 1); // Text color (white)
	indexTextSep->addChild(textColorMat);


	SoText2* indexText = new SoText2();
	SoTranslation* textTranslation = new SoTranslation();

	for (int hv = 0; hv < mesh->samples.size(); ++hv)
	{
		// Create a new SoText2 node for each index
		SoText2* currentIndexText = new SoText2();
		currentIndexText->string = SbString(std::to_string(hv).c_str());

		// Position the text slightly offset from the point
		SoTranslation* currentTextTranslation = new SoTranslation();
		currentTextTranslation->translation.setValue(
			mesh->vertices[mesh->samples[hv]]->coords[0] + 0.025f, // Adjust offset as needed
			mesh->vertices[mesh->samples[hv]]->coords[1] + 0.025f,
			mesh->vertices[mesh->samples[hv]]->coords[2] + 0.025f
		);

		SoSeparator* textTransformSep = new SoSeparator();
		textTransformSep->addChild(currentTextTranslation);
		textTransformSep->addChild(currentIndexText);
		indexTextSep->addChild(textTransformSep);
	}
	res->addChild(indexTextSep);

	res->addChild(highlightedVertexSep);*/

	return res;
}
