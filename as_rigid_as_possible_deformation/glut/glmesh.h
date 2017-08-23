//
//    File: glmesh.h
//
//    (C) 2000-2008 Helmut Cantzler
//
//    Licensed under the terms of the Lesser General Public License.
//
//    Shayan:
//    All the GUI is described in this header file.

#ifndef _OPENGL_FUNCTIONS_H
#define _OPENGL_FUNCTIONS_H

#include <GL/glut.h>
#include <GL/glu.h>
#include <stdio.h>
#include <cassert>

#ifndef WIN32
   #include <unistd.h>
#endif
 
#include "mesh.h"
#include "arap.h"

// Types

enum DisplayType {VERTICES=0, SOLID, WIRE, FRONTLINES, NOTHING};
enum ActionType {NONE, TRANSLATE, ROTATE, ZOOM, SELECT, ARAP_MOVE, ARAP_ROTATE};

struct GLMeshSettings
{
  Mesh *mesh;
  Deformation *deform;

  // Viewpoint data (Trackball, Translation, Zoom, Clipping)
  // ModelView Matrix
  GLfloat tbTransform[4][4];
  // Inverse of ModelView Matrix
  GLfloat tbTransformInv[16];
  // These are present in order to imitate zooming with an orthographic camera.
  GLfloat xShift, yShift, zShift;
  // Far and near clipping.
  GLfloat clipping;
 
  // Viewmode (surface, edges+surface, edges, vertices)
  enum DisplayType meshDisplayMode;
 
  // Light intensity
  GLfloat lightBrightness;
 
  // Background colour
  GLfloat red, green, blue;

  // Flags
  bool displayNormals, displayBoundingBox;
  bool enableCutBackFaces, enableAspectRatio;
  // Show anchor vertices?
  bool showSelection;
};

// Colours
constexpr unsigned int N_COLORS = 8;
enum colors {BLACK=0, BLUE, GREEN, RED, YELLOW, PURPLE, CYAN, WHITE};

const GLfloat ambientLightSource[] = {0.0, 0.0, 0.0, 1.0};

extern GLMeshSettings settings;

extern float surfaceNormalLength;
extern GLint windowWidth;
extern GLint windowHeight;
extern GLint windowX;
extern GLint windowY;
extern int lastMouseX, lastMouseY;
extern int selectMouseX, selectMouseY;
extern enum ActionType actionMode;
extern GLfloat tbLastPosition[3];
extern GLfloat tbAxis[3];
extern int activeGroup;
extern int translationZoom;
extern int selectMoveRot;
extern int selectBack;

void draw(void);
void reshape(int width, int height);
void setTitle(void);
void key(unsigned char k, int x, int y);
void arrowKey(int k, int x, int y);
void mouseButton(int button, int state, int x, int y);
void mouseDrag(int x, int y);
void mainMenu(int value);
void displayMenu(int value);
void openglInit(GLMeshSettings s);
void openglStart(void);

void checkForTraps(const char *string, int value);
void setMaterialColor(GLenum face, GLenum pname, int color);
void setMaterialColorGrey(GLenum face, GLenum pname, int b);
void displayNormals(void);
void displayBoundingBox(void);
void displayPoints(list<Vertex*> *vertices);
void displayMesh(void);
void displayMeshWithStencil(void);
void tbPointToVector(int x, int y, float v[3]);
void initLighting(void);

int width(void);
int height(void);

void displaySelectionBox(void);
void gluScrToFrust(int, int, float*, float*);
void selectVertices(int x0, int y0, int x1, int y1, int iselect=0);
void gluUnproject(int windx, int windy, float ans[3]);
template<typename T> void gluMatMultVec(const T mm[16], const T v[4], T ans[4]);
template<typename T> bool gluInvertMatrix(const T m[16], T invOut[16]);

#endif
