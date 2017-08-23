//
//    File: main.cc
//
//    (C) 2000-2008 Helmut Cantzler
//
//    Licensed under the terms of the Lesser General Public License.
//

#include "mesh_io.h"
#include "arap.h"

#include "glmesh.h"

void setNameAndPath(char *n, Mesh *mesh)
{
	char *lastSlash;

#ifndef WIN32
	lastSlash = strrchr(n, '/');
#else
	lastSlash = strrchr(n, '\\');
#endif

	if (lastSlash == 0)
	{
		mesh->setName(n);
		mesh->setPath("");
	}
	else
	{
		lastSlash++;
		mesh->setName(lastSlash);
		*lastSlash = '\0';
		mesh->setPath(n);
	}
}

void initSettings(GLMeshSettings &settings)
{
	settings.mesh = NULL;
	settings.deform = NULL;
	settings.meshDisplayMode = SOLID;
	settings.xShift = settings.yShift = settings.zShift = 0.0;
	settings.clipping = 0.0;
	settings.lightBrightness = 1.0;
	settings.red = settings.green = settings.blue = 0.0;
	settings.displayNormals = false;
	settings.displayBoundingBox = false;
	settings.enableCutBackFaces = false;
	settings.enableAspectRatio = true;
	settings.showSelection = true;

	for (int i=0; i < 4; i++)
		for (int j=0; j < 4; j++)
		{
			settings.tbTransformInv[i+4*j]= i == j ? 1.0 : 0.0;
			settings.tbTransform[i][j]= i == j ? 1.0 : 0.0;
		}
}

int main(int argc, char **argv)
{
	PetscErrorCode ierr;

	glutInit(&argc, argv);
	ierr=PetscInitialize(&argc,&argv, NULL, NULL); CHKERRQ(ierr);

	printf("\nMesh Viewer 0.3.3\n\n");

	if (argc < 2)
	{
		printf("Syntax: %s <mesh file> \n\n", argv[0]);
		exit(1);
	}

	/*
	 * Create the mesh and settings.
	 */
	GLMeshSettings settings;
	initSettings(settings);
	settings.mesh = new Mesh;

	/*
	 * Read the mesh
	 */
	MeshIO io(settings.mesh);
	io.read_auto(argv[1]);
	setNameAndPath(argv[1], settings.mesh);

	printf("[hooshi] mesh name and path: %s, %s", settings.mesh->getName(), settings.mesh->getPath());
	printf("Mesh: %d vertices, %d triangles, and %d edges \n",
			settings.mesh->numberOfVertices(),settings.mesh->numberOfTriangles(),settings.mesh->numberOfEdges() );


	/*
	 * Map vertices and triangle centroids to mesh centroid
	 */
	printf("Map mesh to centre ...\n");
	settings.mesh->moveToCentre();
	printf("Scale mesh ...\n");
	settings.mesh->scaleIntoNormalSphere();


	/*
	 * Create the deformation
	 */
	settings.mesh->createEdges();
	settings.deform = new Deformation(settings.mesh); CHKERRQ(settings.deform->ierr);

	openglInit(settings);

	printf("Start displaying ...\n");

	openglStart();

	// The program will never reach here.
	// Control things from key(ESCAPE)

	return 0;
}
