//
//    File: glmesh.cc
//
//    (C) 2000-2008 Helmut Cantzler
//
//    Licensed under the terms of the Lesser General Public License.
//

#include "glmesh.h"
#include "petscsys.h"

// Texts
#define ASCII_ESCAPE 27

struct MenuItem
{
	const char* text;
	const char key;
};

std::vector<MenuItem>   mainMenuText =
{
		{"Toggled selecting back facing vertices (a)", 'a'},
		{"Toggled viewing anchor vertices (s)", 's'},
		{"Toggled normals (n)", 'n'},
		{"Toggled bounding box (b)", 'b'},
		{"Lighting up (u)", 'u'},
		{"Lighting down (d)", 'd'},
		{"Clipping up (+)", '+'},
		{"Clipping down (-)", '-'},
		{"Full screen (F)", 'F'},
		{"Quit (esc)", ASCII_ESCAPE}
};

std::vector<MenuItem>  displayMenuText =
{
		{"Solid (1)", '1'},
		{"Frontlines (2)", '2'},
		{"Wireframe (3)", '3'},
		{"Points (4)", '4'}
};

std::vector<MenuItem> actionMenuText =
{
		{"Select the anchors (,)", ','},
		{"Move the anchors (.)", '.'},
		{"Deselect All Anchor Group ({)", '{'},
		{"Deselect This Anchor Group (})", '}'}
};

char title[100];

GLMeshSettings settings;

float surfaceNormalLength;
GLint windowWidth = 500;
GLint windowHeight = 500;
GLint windowX = 100;
GLint windowY = 100;
int lastMouseX, lastMouseY;
int selectMouseX, selectMouseY;
enum ActionType actionMode = NONE;
GLfloat tbLastPosition[3];
GLfloat tbAxis[3];
int translationZoom = 0;
int activeGroup = -1;
int selectBack = 0;
int selectMoveRot = 0;

////////////////// GLUT handlers ////////////////// 

void draw(void)
{
	// Clear the current display
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );


	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-2.1 - settings.xShift - settings.zShift,
			2.1 - settings.xShift + settings.zShift,
			-2.1 - settings.yShift - settings.zShift,
			2.1 - settings.yShift + settings.zShift,
			-3.1 - settings.zShift - settings.clipping, 3.1);


	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMultMatrixf((GLfloat*) settings.tbTransform);

	// display main object (typically a mesh)
	switch (settings.meshDisplayMode)
	{
	case VERTICES:
		displayPoints(settings.mesh->getVertices());
		break;


	case SOLID:
		// Set polygon mode
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		displayMesh();
		break;

	case FRONTLINES:
		displayMeshWithStencil();
		break;

	case WIRE:
		// Set polygon mode
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		displayMesh();
		break;

	default:
		printf("Die a horrible death here! \n");
		assert(0); throw;
	}

	checkForTraps("Draw Polygons", 1);


	// Draw the surface normals
	if (settings.displayNormals)
		displayNormals();

	checkForTraps("Draw Surface Normals",1);

	// Draw the bounding box
	if (settings.displayBoundingBox)
		displayBoundingBox();

	checkForTraps("Draw Bounding Box", 1);

	// Draw selection Box
	if ((actionMode == SELECT) )
		displaySelectionBox();

	checkForTraps("Draw Selection Box", 1);

	glMatrixMode(GL_MODELVIEW);

	glutSwapBuffers(); // Needed when double buffering
}

void gluScrToFrust(int xi, int yi, float* xo, float* yo)
{
	int max = ( ( width() > height() ) ? width() : height() );

	if (settings.enableAspectRatio)
	{
		*xo = (2.0 * xi - width()  ) / float(max);
		*yo = (height() - 2.0 * yi  ) / float(max);
	}
	else
	{
		*xo = (2.0 * xi - width() ) / float(width());
		*yo = (height() - 2.0 * yi) / float(height());
	}
}

void reshape(int width, int height)
{
	int max = width > height ? width : height;

	if (settings.enableAspectRatio)
		glViewport((width-max)/2, (height-max)/2, max, max);
	else
		glViewport(0, 0, width, height);

	windowWidth = width;
	windowHeight = height;
}

void setTitle(void)
{
	char str[400];

	strcpy(title, "Mesh Viewer -- ");

	// Active Group or rotation
	if (activeGroup == -1)
		strcat(title, "Rotate");
	else if (!selectMoveRot)
	{
		sprintf(str, "Anchor verts(%d) -- Select(%s)", activeGroup+1, (selectBack ? "b" : "nb") );
		strcat(title, str);
	}
	else if (selectMoveRot == 1)
	{
		sprintf(str, "Anchor verts(%d) -- Move", activeGroup+1 );
		strcat(title, str);
	}
	else
	{
		sprintf(str, "Anchor verts(%d) -- Rotate", activeGroup+1 );
		strcat(title, str);
	}

	strcat(title, " -- ");

	if (translationZoom)
		strcat(title, "Translation");
	else
		strcat(title, "Zoom");
	strcat(title, " -- ");

	switch (settings.meshDisplayMode)
	{
	case SOLID:
		strcat(title, "Solid");
		break;
	case FRONTLINES:
		strcat(title, "Frontlines");
		break;
	case WIRE:
		strcat(title, "Wireframe");
		break;
	case VERTICES:
		strcat(title, "Points");
		break;
	default:
		printf("die a horrible death here! \n");
		assert(0); throw;
	}

	glutSetWindowTitle(title);
}


void key(unsigned char k, int x, int y)
{
	PetscErrorCode ierr;

	switch (k)
	{
	// Quit --------------------------------------
	case ASCII_ESCAPE: // Quit (Esc)
		if(settings.deform)
		{
			ierr=settings.deform->destroy();CHKERRV(ierr);
			delete settings.deform;
		}
		if(settings.mesh) delete settings.mesh;
		ierr=PetscFinalize();CHKERRV(ierr);
		exit(0);
		break;
	// Display Modes-------------------------------

	case 'q':
		settings.meshDisplayMode=SOLID;
		setTitle();
		break;
	case 'w':
		settings.meshDisplayMode=FRONTLINES;
		setTitle();
		break;
	case 'e':
		settings.meshDisplayMode=WIRE;
		setTitle();
		break;
	case 'r':
		settings.meshDisplayMode=VERTICES;
		setTitle();
		break;

	// Selection Groups ---------------------------
	case '1':
	case '2':
	case '3':
	case '4':
	case '5':
	case '6':
	case '7':
		activeGroup = k - '1';
		setTitle();
		break;

	case '0':
		activeGroup = -1;
		setTitle();
		break;

	// Actions -------------------------------------

	case '/': // Rotate
	{
		selectMoveRot = 2;
		setTitle();
		break;
	}

	case '.': // Move
	{
		selectMoveRot = 1;
		setTitle();
		break;
	}

	case ',': // Select
	{
		selectMoveRot = 0;
		setTitle();
		break;
	}

	case 't': // Toggled translation/zoom (t)
	{
		translationZoom=!translationZoom;
		setTitle();
		break;
	}

	case 'a': // Toggle back selection
		selectBack = !selectBack;
		setTitle();
		break;

	case 's': // Show selection
		settings.showSelection = !settings.showSelection;
		glutPostRedisplay();
		break;

	case '{': // Deselect all
	{
		for (int i = 0 ; i < settings.mesh->n_select() ; i++) settings.mesh->clearSelection(i);
		glutPostRedisplay();
		break;
	}

	case '}': // Deselect This
	{
		if (activeGroup >= 0)
		{
			settings.mesh->clearSelection(activeGroup);
			glutPostRedisplay();
		}
		break;
	}

	case '\\':
		settings.deform->resetMeshPositions();
		break;

	// Other guys  --------------------------------
	case 'c':
		settings.enableCutBackFaces=!settings.enableCutBackFaces;
		if (settings.enableCutBackFaces)
			glEnable(GL_CULL_FACE);
		else
			glDisable(GL_CULL_FACE);
		break;


	case 'F': // Full screen (f)
		glutFullScreen();
		break;

	case 'n': // Toggled normals (n)
		settings.displayNormals=!settings.displayNormals;
		break;

	case 'b': // Toggled bounding box (b)
		settings.displayBoundingBox=!settings.displayBoundingBox;
		break;

	case '=': case '+': // Clipping up (+)
		settings.clipping+=0.1;
		break;

	case '-': // Clipping down (-)
		settings.clipping-=0.1;
		break;

	case 'u': // Light up (u)
	{
		settings.lightBrightness+=settings.lightBrightness*0.1;
		GLfloat source[] = {settings.lightBrightness,
				settings.lightBrightness,
				settings.lightBrightness,
				settings.lightBrightness};

		glLightfv(GL_LIGHT0, GL_DIFFUSE, source);
		glLightfv(GL_LIGHT1, GL_DIFFUSE, source);
		glLightfv(GL_LIGHT2, GL_DIFFUSE, source);
		glLightfv(GL_LIGHT3, GL_DIFFUSE, source);
	}
	break;

	case 'd': // Light down (d)
	{
		settings.lightBrightness-=settings.lightBrightness/1.1*0.1;
		GLfloat source[] = {settings.lightBrightness,
				settings.lightBrightness,
				settings.lightBrightness,
				settings.lightBrightness};

		glLightfv(GL_LIGHT0, GL_DIFFUSE, source);
		glLightfv(GL_LIGHT1, GL_DIFFUSE, source);
		glLightfv(GL_LIGHT2, GL_DIFFUSE, source);
		glLightfv(GL_LIGHT3, GL_DIFFUSE, source);

		break;
	}
	} // End of switch

	glutPostRedisplay();
}



void mouseButton(int button, int state, int x, int y)
{
	switch(state)
	{
	case GLUT_DOWN:
	{
		lastMouseX = x;
		lastMouseY = y;

		switch(button)
		{

		case GLUT_LEFT_BUTTON:
		{
			if(activeGroup == -1)
			{
				actionMode = ROTATE;
				tbPointToVector(x, y, tbLastPosition);
			}
			else if (selectMoveRot == 0)
			{
				actionMode = SELECT;
			}
			else if (selectMoveRot == 1)
			{
				actionMode = ARAP_MOVE;
			}
			else if (selectMoveRot == 2)
			{
				actionMode = ARAP_ROTATE;
				tbPointToVector(x, y, tbLastPosition);
			}
			else
			{
				printf("Die a horrible death! \n");
				throw;
			}
			break;
		}

		case GLUT_MIDDLE_BUTTON:
		{
			if (translationZoom)
				actionMode = TRANSLATE;
			else
				actionMode = ZOOM;
			break;
		}

		} // End of switch button

		break;
	} // end of case GLUT_DOWN

	case GLUT_UP:
	{
		if(actionMode == SELECT)
		{
			assert(activeGroup >= 0);
			selectVertices(lastMouseX, lastMouseY, selectMouseX, selectMouseY, activeGroup);
		}
		actionMode = NONE;
		glutPostRedisplay();
		break;
	} // end of case GLUT_UP

	} // end of switch(state)
}

void mouseDrag(int x, int y)
{
	switch (actionMode)
	{
	case ZOOM:
	{
		settings.zShift += (y-lastMouseY)/(double)windowHeight*3;
		lastMouseX = x;
		lastMouseY = y;
		glutPostRedisplay();
		break;
	} // End of case ZOOM

	case TRANSLATE:
	{
		settings.xShift += (x-lastMouseX)/(double)windowWidth*3;
		settings.yShift -= (y-lastMouseY)/(double)windowHeight*3;
		lastMouseX = x;
		lastMouseY = y;
		glutPostRedisplay();
		break;
	} // End of case ZOOM

	case ROTATE:
	{
		float currentPosition[3], angle, dx, dy, dz;

		tbPointToVector(x, y, currentPosition);

		dx = currentPosition[0] - tbLastPosition[0];
		dy = currentPosition[1] - tbLastPosition[1];
		dz = currentPosition[2] - tbLastPosition[2];
		angle = 90.0 * sqrt(dx * dx + dy * dy + dz * dz) * 2;

		// tbLastPosition \cross tbNewPosition
		tbAxis[0] = tbLastPosition[1] * currentPosition[2] -
				tbLastPosition[2] * currentPosition[1];
		tbAxis[1] = tbLastPosition[2] * currentPosition[0] -
				tbLastPosition[0] * currentPosition[2];
		tbAxis[2] = tbLastPosition[0] * currentPosition[1] -
				tbLastPosition[1] * currentPosition[0];

		tbLastPosition[0] = currentPosition[0];
		tbLastPosition[1] = currentPosition[1];
		tbLastPosition[2] = currentPosition[2];

		glPushMatrix();
		glLoadIdentity();
		glRotatef(angle, tbAxis[0], tbAxis[1], tbAxis[2]);
		glMultMatrixf((GLfloat*) settings.tbTransform);
		glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat *)settings.tbTransform);
		gluInvertMatrix((GLfloat *)settings.tbTransform, settings.tbTransformInv);
		glPopMatrix();

		glutPostRedisplay();
		break;
	} // End of case ROTATE

	case SELECT:
	{
		selectMouseX = x;
		selectMouseY = y;
		glutPostRedisplay();
		break;
	}

	case ARAP_MOVE:
	{
		float pr[3], cr[3];

		gluUnproject(lastMouseX, lastMouseY, pr);
		gluUnproject(x, y, cr);
		MathVector delta(cr[0] - pr[0], cr[1] - pr[1], cr[2] - pr[2]);
		if(delta.length() < 0.2) return;

		vector<int> anchor;
		vector<float> newpos;

		/*
		 * Move all the selected points.
		 */
		for (int i = 0 ; i < settings.mesh->n_select() ; i++)
		{
			set<Vertex *> *setv = settings.mesh->getSelectedVertices(i);
			Vertex *v;

			if ( i == activeGroup)
			{
				for (auto iv = setv->begin() ; iv != setv->end() ; ++iv)
				{
					v = *iv;
					anchor.push_back(v->name);
					newpos.push_back(v->x() + delta.v[0]);
					newpos.push_back(v->y() + delta.v[1]);
					newpos.push_back(v->z() + delta.v[2]);
				}
			} // end of active selection
			else
			{
				for (auto iv = setv->begin() ; iv != setv->end() ; ++iv)
				{
					v = *iv;
					anchor.push_back(v->name);
					newpos.push_back(v->x());
					newpos.push_back(v->y());
					newpos.push_back(v->z());
				}
			} // end of non active selection
		} // End off loop over selected verts

		/*
		 * Deform
		 */
		settings.deform->deform(anchor, newpos);

		/*
		 * Update mouse
		 */
		lastMouseX=x;
		lastMouseY=y;

		/*
		 * Redisplay
		 */
		glutPostRedisplay();
		break;
	} // End of case ARAP_MOVE

	case ARAP_ROTATE:
	{
		break;
		set<Vertex *> *setv;

		/*
		 * Find the rotation matrix
		 */
		float currentPosition[3], angle, dx, dy, dz;
		GLfloat rotationMatrix[16];

		tbPointToVector(x, y, currentPosition);

		dx = currentPosition[0] - tbLastPosition[0];
		dy = currentPosition[1] - tbLastPosition[1];
		dz = currentPosition[2] - tbLastPosition[2];

		if(sqrt(dx*dx + dy*dy + dz*dz) < 0.2) break;

		angle = 90.0 * sqrt(dx * dx + dy * dy + dz * dz) * 2;

		// tbLastPosition \cross tbNewPosition
		tbAxis[0] = tbLastPosition[1] * currentPosition[2] -
				tbLastPosition[2] * currentPosition[1];
		tbAxis[1] = tbLastPosition[2] * currentPosition[0] -
				tbLastPosition[0] * currentPosition[2];
		tbAxis[2] = tbLastPosition[0] * currentPosition[1] -
				tbLastPosition[1] * currentPosition[0];

		tbLastPosition[0] = currentPosition[0];
		tbLastPosition[1] = currentPosition[1];
		tbLastPosition[2] = currentPosition[2];


		float cent[3] = {0,0,0};
		setv = settings.mesh->getSelectedVertices(activeGroup);
		for (auto iv = setv->begin() ; iv != setv->end() ; ++iv)
		{
			Vertex *v = *iv;
			cent[0] += v->x();
			cent[1] += v->y();
			cent[2] += v->z();
		}
		cent[0] /= setv->size();
		cent[1] /= setv->size();
		cent[2] /= setv->size();

		glPushMatrix();
		glLoadIdentity();
		glMultMatrixf((GLfloat*)settings.tbTransform);
		glTranslatef(cent[0],cent[1],cent[2]);
		glRotatef(angle, tbAxis[0], tbAxis[1], tbAxis[2]);
		glTranslatef(-cent[0],-cent[1],-cent[2]);
		glMultMatrixf((GLfloat*)settings.tbTransformInv);
		glGetFloatv(GL_MODELVIEW_MATRIX, rotationMatrix);
		glPopMatrix();

		/*
		 * Set the anchors
		 */
		vector<int> anchor;
		vector<float> newpos;
		/*
		 * Move all the selected points.
		 */
		for (int i = 0 ; i < settings.mesh->n_select() ; i++)
		{
			setv = settings.mesh->getSelectedVertices(i);
			Vertex *v;

			if ( i == activeGroup)
			{
				for (auto iv = setv->begin() ; iv != setv->end() ; ++iv)
				{
					v = *iv;
					anchor.push_back(v->name);
					float newXYZ[4];
					float oldXYZ[4] = {v->x(), v->y(), v->z(), 1. };
					gluMatMultVec(rotationMatrix, oldXYZ, newXYZ);
					newpos.push_back(newXYZ[0]);
					newpos.push_back(newXYZ[1]);
					newpos.push_back(newXYZ[2]);
				}
			} // end of active selection
			else
			{
				for (auto iv = setv->begin() ; iv != setv->end() ; ++iv)
				{
					v = *iv;
					anchor.push_back(v->name);
					newpos.push_back(v->x());
					newpos.push_back(v->y());
					newpos.push_back(v->z());
				}
			} // end of non active selection
		} // End off loop over selected verts

		/*
		 * Deform
		 */
		settings.deform->deform(anchor, newpos);

		/*
		 * Update mouse
		 */
		lastMouseX=x;
		lastMouseY=y;

		/*
		 * Redisplay
		 */
		glutPostRedisplay();
		break;
	} // End of case (ARAP_ROTATE)

	default: ;
	// do nothing
	}
}

void arrowKey(int k, int x, int y)
{
	if(activeGroup < 0) return;

	PetscErrorCode ierr;
	set<Vertex *> *setv;
	GLfloat rotationMatrix[16];
	const float angle = 5.;

	/*
	 * Find the rotation matrix
	 */

	// FInd center
	float cent[3] = {0,0,0};
	setv = settings.mesh->getSelectedVertices(activeGroup);
	for (auto iv = setv->begin() ; iv != setv->end() ; ++iv)
	{
		Vertex *v = *iv;
		cent[0] += v->x();
		cent[1] += v->y();
		cent[2] += v->z();
	}
	cent[0] /= setv->size();
	cent[1] /= setv->size();
	cent[2] /= setv->size();

	// create the matrix
	glPushMatrix();
	glLoadIdentity();
	glMultMatrixf((GLfloat*)settings.tbTransform);
	glTranslatef(cent[0],cent[1],cent[2]);
	switch (k)
	{
	case GLUT_KEY_UP:
		glRotatef(angle, 0+1, 0, 0 );
		break;
	case GLUT_KEY_DOWN:
		glRotatef(angle, 0-1,0 , 0);
		break;
	case GLUT_KEY_LEFT:
		glRotatef(angle, 0, 0-1, 0);
		break;
	case GLUT_KEY_RIGHT:
		glRotatef(angle, 0,0+1 ,0 );
		break;
	case GLUT_KEY_PAGE_UP:
		glRotatef(angle, 0, 0, 0-1);
		break;
	case GLUT_KEY_PAGE_DOWN:
		glRotatef(angle, 0, 0, 0+1);
		break;
	}
	glTranslatef(-cent[0],-cent[1],-cent[2]);
	glMultMatrixf((GLfloat*)settings.tbTransformInv);
	glGetFloatv(GL_MODELVIEW_MATRIX, rotationMatrix);
	glPopMatrix();


	/*
	 * Set the anchors
	 */
	vector<int> anchor;
	vector<float> newpos;
	/*
	 * Move all the selected points.
	 */
	for (int i = 0 ; i < settings.mesh->n_select() ; i++)
	{
		setv = settings.mesh->getSelectedVertices(i);
		Vertex *v;

		if ( i == activeGroup)
		{
			for (auto iv = setv->begin() ; iv != setv->end() ; ++iv)
			{
				v = *iv;
				anchor.push_back(v->name);
				float newXYZ[4];
				float oldXYZ[4] = {v->x(), v->y(), v->z(), 1. };
				gluMatMultVec(rotationMatrix, oldXYZ, newXYZ);
				newpos.push_back(newXYZ[0]);
				newpos.push_back(newXYZ[1]);
				newpos.push_back(newXYZ[2]);
			}
		} // end of active selection
		else
		{
			for (auto iv = setv->begin() ; iv != setv->end() ; ++iv)
			{
				v = *iv;
				anchor.push_back(v->name);
				newpos.push_back(v->x());
				newpos.push_back(v->y());
				newpos.push_back(v->z());
			}
		} // end of non active selection
	} // End off loop over selected verts

	/*
	 * Deform
	 */
	ierr=settings.deform->deform(anchor, newpos);CHKERRV(ierr);

	/*
	 * Redisplay
	 */
	glutPostRedisplay();

}

void mainMenu(int value)
{
	key(mainMenuText[value].key,0,0);
}

void displayMenu(int value)
{
	key(displayMenuText[value].key,0,0);
}

void actionMenu(int value)
{
	key(actionMenuText[value].key,0,0);
}

////////////////// Intialisation ////////////////// 

void openglInit(GLMeshSettings s)
{
	settings=s;

	settings.meshDisplayMode=NOTHING;
	if (settings.mesh != NULL)
	{
		if (settings.mesh->numberOfTriangles() > 0)
			settings.meshDisplayMode=SOLID;
		else if (settings.mesh->numberOfVertices() > 0)
			settings.meshDisplayMode=VERTICES;
		surfaceNormalLength = 0.02 + settings.mesh->averageTriangleSize()/10;
	}
	else
	{
		printf("Cannot proceed with NULL mesh. \n");
		assert(0 && "Cannot proceed with NULL mesh"); throw;
	}

	actionMode=NONE;

	// Initialise the display window
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB | GLUT_STENCIL);
	glutInitWindowSize(windowWidth, windowHeight);
	glutInitWindowPosition(windowX, windowY);
	glutCreateWindow("");
	setTitle();

	// Set background colour and clear
	glClearColor(0.0, 0.0, 0.0, 0.0);  // Black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Initialise light source
	initLighting();

	// Other initialisations
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glLineWidth(1.5);
	glPointSize(1.5);
	//glEnable(GL_POINT_SMOOTH);

	if (settings.enableCutBackFaces)
		glEnable(GL_CULL_FACE);
	else
		glDisable(GL_CULL_FACE);

	// Enable Flat Shading
	glShadeModel(GL_FLAT);
	// Really Nice Perspective Calculations
	//glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

	// Initialise the modelview transformation
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// Make the menu

	glutCreateMenu(mainMenu);
	glutCreateMenu(displayMenu);
	glutCreateMenu(actionMenu);

	glutSetMenu(2);
	for (unsigned int i=0; i < displayMenuText.size(); i++)
		glutAddMenuEntry(displayMenuText[i].text, i);

	glutSetMenu(3);
	for (unsigned int i=0; i < actionMenuText.size(); i++)
		glutAddMenuEntry(actionMenuText[i].text, i);

	glutSetMenu(1);
	glutAddSubMenu("Display", 2);
	glutAddSubMenu("Action", 3);
	for (unsigned int i=0; i < mainMenuText.size(); i++)
		glutAddMenuEntry(mainMenuText[i].text, i);

	glutAttachMenu(GLUT_RIGHT_BUTTON);
}

void openglStart(void)
{
	// Register callback handlers
	glutDisplayFunc(draw);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(key);
	glutMouseFunc(mouseButton);
	glutMotionFunc(mouseDrag);
	glutSpecialFunc(arrowKey);
	glutIdleFunc(0);

	// Main (infinite) loop
	glutMainLoop();
}

int width(void)
{
	return windowWidth;
}

int height(void)
{
	return windowHeight;
}
