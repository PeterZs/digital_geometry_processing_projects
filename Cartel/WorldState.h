/* Copyright (c) Russell Gillette
 * December 2013
 *
 * CPSC Computer Graphics 2014
 * University of British Columbia
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
 * documentation files (the "Software"), to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and
 * to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF
 * OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */


/* == WorldState.h ==
 *
 * global data specific to the program and
 * related GPU operations
 */

#ifndef WORLD_STATE_H
#define WORLD_STATE_H

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <vector>
#include <cstdio>
#include <math.h>

#ifndef GLM_FORCE_RADIANS 
#define GLM_FORCE_RADIANS 
#endif
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "TextureTypes.h"

// Lighting state
struct LightInfo
{
    glm::vec4 LPosition;
    glm::vec3 La;       // ambient light
    glm::vec3 Ld;       // diffuse light
    glm::vec3 Ls;       // specular light
    glm::vec3 Intensity;// light intensity

    // some default initialization parameters
    LightInfo()
    : LPosition(10.0, 6.0, 4.0, 1.0),
    La(0.9, 0.9, 0.9),
    Ld(1.0, 1.0, 1.0),
    Ls(1.0, 1.0, 1.0),
    Intensity(1.0,1.0,1.0){}
};

// Material and Color state
struct MaterialInfo
{
    glm::vec3 Ka;
    glm::vec3 Kd;
    glm::vec3 Ks;
    float Shininess;

    // some default initialization parameters
    MaterialInfo()
        : Ka(0.5, 0.5, 0.5),
          Kd(0.7, 0.7, 0.7),
          Ks(0.9, 0.9, 0.9),
          Shininess(100.0f) {}
    void custom_color(glm::vec3 color) {Ka = color;}
};

/**
 * Encapsulation of common global variables as well as related GPU operations
 */
class WorldState
{
public:
    WorldState();
    ~WorldState();

    // store and set current active shaders
    void useProgram(unsigned int i);
    void setProgram(unsigned int i, GLuint program);
    GLuint getCurrentProgram() {return shaders[currentProgram];}

    // send transform matrices, lights, materials to the GPU
    void loadTransforms();                          // load global matrices
    void loadTransforms(glm::mat4 M);              // load custom MV matrix, use global P
    void loadTransforms(glm::mat4 M, glm::mat4 V, glm::mat4 P); 
    void loadObjectTransforms(glm::mat4 oM);       // load the transformations of an object (multiply by global)
    void loadLights();                              // load all lights
    void loadLight(unsigned int i);                 // load specific light
    void loadMaterials();                           // load all materials
    void loadMaterial(unsigned int i);              // load specific material
    void loadColorMaterial(glm::vec4 color);        // load the passed in color as material0
    void loadTexture(unsigned int i);
    void loadTextures();

    void resetProjection(float aspectRatio);
    void resetView();
    void resetModel();

    void updateView(float vTheta, float vPhi, float vDepth, glm::vec3 vPan);

    //default initialized to the identity matrix
    glm::mat4 model;
    glm::mat4 view;
    glm::mat4 projection;
    glm::vec3 cameraOrigin;


    std::vector<LightInfo>    lights;
    std::vector<MaterialInfo> materials;
    std::vector<GLuint>       shaders;
    std::vector<Texture>      textures;

private:
    unsigned int currentProgram;

};

#endif // WORLD_STATE_H