//
//    File: glmesh_common.cc
//
//    (C) 2000-2008 Helmut Cantzler
//
//    Licensed under the terms of the Lesser General Public License.
//

#include "glmesh.h"

#include <iostream>
#include <iomanip>

void checkForTraps(const char *string, int value)
{
	GLenum errorCode = glGetError();
	if( errorCode!= GL_NO_ERROR)
	{
		printf("\n\n\n");
		printf("During processing %s (%d) an internal error occurred!\n",
				string, value);
		printf("OpenGL error: %d ", errorCode);
		printf("(%s)\n", gluErrorString(errorCode));
		exit(0);
	}
}

void setMaterialColor(GLenum face, GLenum pname, int color)
{
	static float colors[8][4] = {{0.0, 0.0, 0.0, 0.0},  // black
			{0.0, 0.0, 1.0, 1.0},  // blue
			{0.0, 1.0, 0.0, 1.0},  // green
			{1.0, 0.0, 0.0, 1.0},  // red
			{1.0, 1.0, 0.0, 1.0},  // yellow
			{1.0, 0.0, 1.0, 1.0},  // purple
			{0.0, 1.0, 1.0, 1.0},  // cyan
			{1.0, 1.0, 1.0, 1.0}}; // white
	float b;
	int c, j;

	color--;
	// b = 0.0 until 0.9
	b = color / 7 < 10 ? color / 7 * 0.1 : 0.9;
	// c = 0 until 7
	c = color != -1 ? color % 7 + 1 : 0;

	for (j=0; j < 3; j++)
		colors[c][j]= colors[c][j] ? 1.0 - b : 0.0;

	glMaterialfv(face, pname, colors[c]);
}

void setMaterialColorGrey(GLenum face, GLenum pname, int b)
{
	static float grey[4] = {1.0, 1.0, 1.0, 1.0};
	int j;

	for (j=0; j < 3; j++)
		grey[j] = 1.0 - 0.15 * b;

	glMaterialfv(face, pname, grey);
}

void displayNormals(void)
{
	list<Triangle*>::iterator it;
	list<Triangle*> *triangles;

	setMaterialColor(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, BLUE);
	glBegin(GL_LINES);

	triangles=settings.mesh->getTriangles();

	if (triangles != NULL)
	{
		for (it=triangles->begin(); it != triangles->end(); it++)
		{				const float *dnorm = (*it)->floatNormal();
		const float *dcent = (*it)->centroid();

		glVertex3fv((GLfloat*) dcent);
		glVertex3f( dcent[0] + surfaceNormalLength*dnorm[0],
				dcent[1] + surfaceNormalLength*dnorm[1],
				dcent[2] + surfaceNormalLength*dnorm[2]);
		}
	}

	glEnd();
}

static void displayBox(float xMin, float yMin, float zMin, float xMax, float yMax, float zMax)
{
	int size;
	glGetIntegerv(GL_LINE_WIDTH, &size);
	glLineWidth(2);

	setMaterialColor(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, WHITE);

	glBegin(GL_LINES);
	glVertex3f(xMin, yMin, zMin);  glVertex3f(xMin, yMin, zMax);
	glVertex3f(xMin, yMin, zMin);  glVertex3f(xMin, yMax, zMin);
	glVertex3f(xMin, yMin, zMin);  glVertex3f(xMax, yMin, zMin);
	glVertex3f(xMin, yMax, zMax);  glVertex3f(xMax, yMax, zMax);
	glVertex3f(xMin, yMax, zMax);  glVertex3f(xMin, yMin, zMax);
	glVertex3f(xMin, yMax, zMax);  glVertex3f(xMin, yMax, zMin);
	glVertex3f(xMax, yMin, zMax);  glVertex3f(xMax, yMax, zMax);
	glVertex3f(xMax, yMin, zMax);  glVertex3f(xMax, yMin, zMin);
	glVertex3f(xMax, yMin, zMax);  glVertex3f(xMin, yMin, zMax);
	glVertex3f(xMax, yMax, zMin);  glVertex3f(xMax, yMax, zMax);
	glVertex3f(xMax, yMax, zMin);  glVertex3f(xMax, yMin, zMin);
	glVertex3f(xMax, yMax, zMin);  glVertex3f(xMin, yMax, zMin);
	glEnd();

	glLineWidth(size);

}
void displayBoundingBox(void)
{
	if (settings.mesh == NULL)
		return;

	const float xMin = settings.mesh->getXMin() - 0.1;
	const float xMax = settings.mesh->getXMax() + 0.1;
	const float yMin = settings.mesh->getYMin() - 0.1;
	const float yMax = settings.mesh->getYMax() + 0.1;
	const float zMin = settings.mesh->getZMin() - 0.1;
	const float zMax = settings.mesh->getZMax() + 0.1;
	displayBox(xMin, yMin, zMin, xMax, yMax, zMax);
}

void displaySelectedPoints()
{


	// Using spheres
	GLUquadric *quad = gluNewQuadric();

	for (int i = 0 ; i < settings.mesh->n_select() ; i++)
	{
		float rad;

		colors co = colors(i+1);
		if (co == CYAN) co = BLACK;
		setMaterialColor(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, co);
		rad = 0.025;

		set<Vertex*> *verts = settings.mesh->getSelectedVertices(i);

		for (auto it=verts->begin(); it != verts->end(); it++)
		{

			Vertex *v = *it;

			glPushMatrix();
			glTranslatef(v->x(), v->y(), v->z());

			gluSphere(quad, rad, 20, 20);
			glPopMatrix();
		}

		/*
		if(verts->size())
		{
			glPushMatrix();
			glTranslatef(settings.mesh->getSelectedCenter(i)[0], settings.mesh->getSelectedCenter(i)[1], settings.mesh->getSelectedCenter(i)[2]);
			glScaled(0.05,0.05,0.05);
			glutSolidIcosahedron();
			glPopMatrix();
		}
		*/
	} // End of for over selection sets

	gluDeleteQuadric(quad);

	// Using points
	/*
	{
		int size;

		set<Vertex*> *verts = settings.mesh->getSelectedVertices(0);

		setMaterialColor(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, RED);
		glGetIntegerv(GL_POINT_SIZE, &size);
		glPointSize(5);

		glBegin(GL_POINTS);
		for (auto iv=verts->begin(); iv != verts->end(); iv++)
		{
			if ( (*iv)->floatNormal() != NULL )
			{
				glNormal3fv( (GLfloat*) (*iv)->floatNormal() );
			}
			glVertex3fv( (GLfloat*) (*iv)->floatData() );
		}
		glEnd();

		glPointSize(size);
	}
	*/

	// Normal for selected vertices
	/*
	{

		const double len = surfaceNormalLength * 2;
		int size;
		set<Vertex*> *verts = settings.mesh->getSelectedVertices(1);

		glGetIntegerv(GL_LINE_WIDTH, &size);
		setMaterialColor(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, GREEN);
		glLineWidth(3);

		glBegin(GL_LINES);
		for (auto it=verts->begin(); it != verts->end(); it++)
		{
			const float *dnorm = (*it)->floatNormal();
			const float *dcent = (*it)->floatData();

			glVertex3fv((GLfloat*) dcent);
			glVertex3f( dcent[0] + len*dnorm[0],
					dcent[1] + len*dnorm[1],
					dcent[2] + len*dnorm[2]);
		}
		glEnd();

		glLineWidth(size);
	}
	*/
}
void displayPoints(list<Vertex*> *vertices)
{
	setMaterialColor(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, YELLOW);

	for (auto iv=vertices->begin(); iv != vertices->end(); iv++)
	{
		//Vertex *v = *iv;

		//glLoadName((*iv)->name);
		glBegin(GL_POINTS);
		if ( (*iv)->floatNormal() != NULL )
		{
			glNormal3fv( (GLfloat*) (*iv)->floatNormal() );
		}
		glVertex3fv( (GLfloat*) (*iv)->floatData() );
		glEnd();
	}

	// Draw the selected vertices
	if(settings.showSelection) displaySelectedPoints();
}


void displayMesh(void)
{
	// Draw the mesh
	list<Triangle*> *triangles = settings.mesh->getTriangles();
	if (triangles == NULL) return;

	setMaterialColor(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, CYAN);

	for (auto it=triangles->begin(); it != triangles->end(); it++)
	{
		glLoadName((*it)->name);
		glBegin(GL_TRIANGLES);
		glNormal3fv( (GLfloat*) (*it)->floatNormal() );
		glVertex3fv( (GLfloat*) (*it)->vertices[0]->floatData() );
		glVertex3fv( (GLfloat*) (*it)->vertices[1]->floatData() );
		glVertex3fv( (GLfloat*) (*it)->vertices[2]->floatData() );
		glEnd();
	}

	// Draw the selected vertices
	if(settings.showSelection) displaySelectedPoints();
}

void displayMeshWithStencil(void)
{
	list<Triangle*>::iterator it;
	list<Triangle*> *triangles;
	int nr;

	glEnable(GL_STENCIL_TEST);
	glClear(GL_STENCIL_BUFFER_BIT);

	glStencilFunc(GL_ALWAYS, 0, 1);
	glStencilOp(GL_INVERT, GL_INVERT, GL_INVERT);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	triangles=settings.mesh->getTriangles();

	if (triangles != NULL)
	{
		nr = CYAN; // BLUE;

		setMaterialColor(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, nr);

		for (it=triangles->begin(); it != triangles->end(); it++)
		{
			glLoadName((*it)->name);
			glBegin(GL_TRIANGLES);
			glVertex3fv( (GLfloat*) (*it)->vertices[0]->floatData() );
			glVertex3fv( (GLfloat*) (*it)->vertices[1]->floatData() );
			glVertex3fv( (GLfloat*) (*it)->vertices[2]->floatData() );
			glEnd();

			glStencilFunc(GL_EQUAL, 0, 1);
			glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			setMaterialColor(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, nr);
			glBegin(GL_TRIANGLES);
			glNormal3fv( (GLfloat*) (*it)->floatNormal() );
			glVertex3fv( (GLfloat*) (*it)->vertices[0]->floatData() );
			glVertex3fv( (GLfloat*) (*it)->vertices[1]->floatData() );
			glVertex3fv( (GLfloat*) (*it)->vertices[2]->floatData() );
			glEnd();

			glStencilFunc(GL_ALWAYS, 0, 1);
			glStencilOp(GL_INVERT, GL_INVERT, GL_INVERT);
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			setMaterialColor(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, RED);
			glBegin(GL_TRIANGLES);
			glNormal3fv( (GLfloat*) (*it)->floatNormal() );
			glVertex3fv( (GLfloat*) (*it)->vertices[0]->floatData() );
			glVertex3fv( (GLfloat*) (*it)->vertices[1]->floatData() );
			glVertex3fv( (GLfloat*) (*it)->vertices[2]->floatData() );
			glEnd();
		}
	}

	glDisable(GL_STENCIL_TEST);

	if(settings.showSelection) displaySelectedPoints();
}

void tbPointToVector(int x, int y, float v[3])
{
	// project x, y onto a hemi-sphere centered within width, height
	v[0] = (2.0 * x - width() ) / width();
	v[1] = (height() - 2.0 * y) / height();
	const float d2 = v[0] * v[0] + v[1] * v[1];
	const float d = sqrt(d2);
	// v[2] = cos( (3.14159265 / 2.0) * ( (d < 1.0) ? d : 1.0) );
	v[2] = ( (d < 1.0) ? 1 - d2 : 0.0 );
	const float a = 1.0 / sqrt(d2 + v[2] * v[2]);
	v[0] *= a;
	v[1] *= a;
	v[2] *= a;
}

void initLighting(void)
{
	// light from a light source
	const GLfloat diffuseLight[] = {0.4, 0.4, 0.4, 1.0};
	// light from no particulat light source
	const GLfloat ambientLight[] = {0.1, 0.1, 0.1, 1.0};
	// light positions for 4 lights
	const GLfloat lightPositions[4][4] = {{ 1.0,  1.0,  0.0, 0.0},
			{-1.0, -1.0,  0.0, 0.0},
			{-0.1, -0.1,  1.0, 0.0},
			{ 0.1,  0.1, -1.0, 0.0}};

	glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);

	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
	glLightfv(GL_LIGHT0, GL_POSITION, lightPositions[0]);
	glEnable(GL_LIGHT0);

	glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuseLight);
	glLightfv(GL_LIGHT1, GL_POSITION, lightPositions[1]);
	glEnable(GL_LIGHT1);

//	glLightfv(GL_LIGHT2, GL_DIFFUSE, diffuseLight);
//	glLightfv(GL_LIGHT2, GL_POSITION, lightPositions[2]);
//	glEnable(GL_LIGHT2);

//	glLightfv(GL_LIGHT3, GL_DIFFUSE, diffuseLight);
//	glLightfv(GL_LIGHT3, GL_POSITION, lightPositions[3]);
//	glEnable(GL_LIGHT3);

	glEnable(GL_LIGHTING);
}

void displaySelectionBox(void)
{

	int size;
	float x0, y0, x1, y1;

	gluScrToFrust(lastMouseX, lastMouseY, &x0, &y0);
	gluScrToFrust(selectMouseX, selectMouseY, &x1, &y1);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glGetIntegerv(GL_LINE_WIDTH, &size);
	glLineWidth(3);

	setMaterialColor(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, RED);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	glBegin(GL_POLYGON);
	glVertex3f(x0, y0, -1);
	glVertex3f(x1, y0, -1);
	glVertex3f(x1, y1, -1);
	glVertex3f(x0, y1, -1);
	glEnd();

	glLineWidth(size);

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

}

template<typename T>
bool gluInvertMatrix(const T m[16], T invOut[16])
{
    T inv[16], det;
    int i;

    inv[0] = m[5]  * m[10] * m[15] -
             m[5]  * m[11] * m[14] -
             m[9]  * m[6]  * m[15] +
             m[9]  * m[7]  * m[14] +
             m[13] * m[6]  * m[11] -
             m[13] * m[7]  * m[10];

    inv[4] = -m[4]  * m[10] * m[15] +
              m[4]  * m[11] * m[14] +
              m[8]  * m[6]  * m[15] -
              m[8]  * m[7]  * m[14] -
              m[12] * m[6]  * m[11] +
              m[12] * m[7]  * m[10];

    inv[8] = m[4]  * m[9] * m[15] -
             m[4]  * m[11] * m[13] -
             m[8]  * m[5] * m[15] +
             m[8]  * m[7] * m[13] +
             m[12] * m[5] * m[11] -
             m[12] * m[7] * m[9];

    inv[12] = -m[4]  * m[9] * m[14] +
               m[4]  * m[10] * m[13] +
               m[8]  * m[5] * m[14] -
               m[8]  * m[6] * m[13] -
               m[12] * m[5] * m[10] +
               m[12] * m[6] * m[9];

    inv[1] = -m[1]  * m[10] * m[15] +
              m[1]  * m[11] * m[14] +
              m[9]  * m[2] * m[15] -
              m[9]  * m[3] * m[14] -
              m[13] * m[2] * m[11] +
              m[13] * m[3] * m[10];

    inv[5] = m[0]  * m[10] * m[15] -
             m[0]  * m[11] * m[14] -
             m[8]  * m[2] * m[15] +
             m[8]  * m[3] * m[14] +
             m[12] * m[2] * m[11] -
             m[12] * m[3] * m[10];

    inv[9] = -m[0]  * m[9] * m[15] +
              m[0]  * m[11] * m[13] +
              m[8]  * m[1] * m[15] -
              m[8]  * m[3] * m[13] -
              m[12] * m[1] * m[11] +
              m[12] * m[3] * m[9];

    inv[13] = m[0]  * m[9] * m[14] -
              m[0]  * m[10] * m[13] -
              m[8]  * m[1] * m[14] +
              m[8]  * m[2] * m[13] +
              m[12] * m[1] * m[10] -
              m[12] * m[2] * m[9];

    inv[2] = m[1]  * m[6] * m[15] -
             m[1]  * m[7] * m[14] -
             m[5]  * m[2] * m[15] +
             m[5]  * m[3] * m[14] +
             m[13] * m[2] * m[7] -
             m[13] * m[3] * m[6];

    inv[6] = -m[0]  * m[6] * m[15] +
              m[0]  * m[7] * m[14] +
              m[4]  * m[2] * m[15] -
              m[4]  * m[3] * m[14] -
              m[12] * m[2] * m[7] +
              m[12] * m[3] * m[6];

    inv[10] = m[0]  * m[5] * m[15] -
              m[0]  * m[7] * m[13] -
              m[4]  * m[1] * m[15] +
              m[4]  * m[3] * m[13] +
              m[12] * m[1] * m[7] -
              m[12] * m[3] * m[5];

    inv[14] = -m[0]  * m[5] * m[14] +
               m[0]  * m[6] * m[13] +
               m[4]  * m[1] * m[14] -
               m[4]  * m[2] * m[13] -
               m[12] * m[1] * m[6] +
               m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] +
              m[1] * m[7] * m[10] +
              m[5] * m[2] * m[11] -
              m[5] * m[3] * m[10] -
              m[9] * m[2] * m[7] +
              m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] -
             m[0] * m[7] * m[10] -
             m[4] * m[2] * m[11] +
             m[4] * m[3] * m[10] +
             m[8] * m[2] * m[7] -
             m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] +
               m[0] * m[7] * m[9] +
               m[4] * m[1] * m[11] -
               m[4] * m[3] * m[9] -
               m[8] * m[1] * m[7] +
               m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] -
              m[0] * m[6] * m[9] -
              m[4] * m[1] * m[10] +
              m[4] * m[2] * m[9] +
              m[8] * m[1] * m[6] -
              m[8] * m[2] * m[5];

    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

    if (det == 0)
        return false;

    det = 1.0 / det;

    for (i = 0; i < 16; i++)
        invOut[i] = inv[i] * det;

    return true;
}
template bool gluInvertMatrix(const double m[16], double invOut[16]);
template bool gluInvertMatrix(const float m[16], float invOut[16]);


template<typename T>
void gluMatMultVec(const T mm[16], const T v[4], T ans[4])
{
	T (*m)[4] = (T(*)[4])mm;

	for (int i = 0 ; i < 4 ; i++)
	{
		ans[i] = 0;
		for (int j = 0 ; j < 4 ; j++)
			ans[i] += m[j][i] * v[j];
	}
}
template void gluMatMultVec(const double mm[16], const double v[4], double ans[4]);
template void gluMatMultVec(const float mm[16], const float v[4], float ans[4]);


template<typename T>
static void gluMatTranspose(const T mm[16], T aans[16])
{
	T (*m)[4] = (T(*)[4])mm;
	T (*ans)[4] = (T(*)[4])aans;

	for (int i = 0 ; i < 4 ; i++)
		for (int j = 0 ; j < 4 ; j++)
			ans[i][j] = m[j][i];
}

template<typename T>
static void gluPrintMat(const T mm[16], const bool spc = true)
{
	T (*m)[4] = (T(*)[4])mm;

	for (int i = 0 ; i < 4 ; i++)
	{
		for (int j = 0 ; j < 4 ; j++)
			cout << setprecision(4) << setw(8) << m[j][i] << "  ";
		cout << endl;
	}
	if (spc) cout << endl;
}


void gluUnproject(int windx, int windy, float ans[3])
{
	// modelview matrix
	GLfloat modelView[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, modelView);

	// projection
	GLfloat projection[16];
	glGetFloatv(GL_PROJECTION_MATRIX, projection);

	// modelviewprojection
	GLfloat MVP[16];
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glMultMatrixf((GLfloat*) projection);
	glMultMatrixf((GLfloat*) modelView);
	glGetFloatv(GL_MODELVIEW_MATRIX, MVP);
	glPopMatrix();

	//Inverse of MVP
	GLfloat INV[16];
	if(!gluInvertMatrix(MVP, INV))
	{
		printf("Could not invert MVP \n");
		return;
	}

	float camvec[4] = {-1, -1, 0, 1};
	float worldvec[4];

	gluScrToFrust(windx, windy, &camvec[0], &camvec[1]);
	gluMatMultVec(INV, camvec, worldvec);

	ans[0] = worldvec[0];
	ans[1] = worldvec[1];
	ans[2] = worldvec[2];
}

void selectVertices(int winx0, int winy0, int winx1, int winy1, int iselect)
{
	// modelview matrix
	GLfloat modelView[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, modelView);

	// projection
	GLfloat projection[16];
	glGetFloatv(GL_PROJECTION_MATRIX, projection);

	// modelviewprojection
	GLfloat MVP[16];
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glMultMatrixf(projection);
	glMultMatrixf(modelView);
	glGetFloatv(GL_MODELVIEW_MATRIX, MVP);
	glPopMatrix();

	// normal (in case there is no scaling should be the same as MVP)
	GLfloat PVM[16], normalMatrix[16];
	gluMatTranspose(MVP, PVM);
	if(!gluInvertMatrix(PVM, normalMatrix))
	{
		printf("could not invert PVM \n");
		return;
	}

	float cambb0[2], cambb1[2];

	// Note: winy0, and winy1 point downwards.
	gluScrToFrust(std::min(winx0, winx1), std::max(winy0, winy1), &cambb0[0], &cambb0[1]);
	gluScrToFrust(std::max(winx0, winx1), std::min(winy0, winy1), &cambb1[0], &cambb1[1]);

	list<Vertex*> * vertices = settings.mesh->getVertices();
	for (auto iv = vertices->begin() ; iv != vertices->end() ; ++iv)
	{
		Vertex *v = *iv;
		float vcam[4];
		float normcam[4];
		float vworld[4] = {v->x(),v->y(),v->z(),1};
		float normworld[4] = {v->floatNormal()[0], v->floatNormal()[1], v->floatNormal()[2], 1};
		gluMatMultVec(MVP, vworld, vcam);
		gluMatMultVec(normalMatrix, normworld, normcam);

		const bool c1 = (vcam[0] > cambb0[0]) && (vcam[1] > cambb0[1]);
		const bool c2 = (vcam[0] < cambb1[0]) && (vcam[1] < cambb1[1]);
		const bool c3 = (selectBack ? true : (normcam[2] < 0) );
		if (c1 && c2 && c3)
		{
			settings.mesh->selectVertex(v, iselect);
		}
	}
}
