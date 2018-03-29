/* Copyright (c) Russell Gillette
 * December 2013
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated
 * documentation files (the "Software"), to deal in the Software without
 * restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and
 * to permit persons to whom the Software is furnished to do so, subject to the
 * following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
 * NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF
 * OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <stdio.h>

#include "WorldState.h"

WorldState::WorldState() : cameraOrigin(0.0f, 2.5f, 10.0f)
{
  // make sure that the world can hold at least 5 of each object
  // before it needs to reallocate memory
  shaders.reserve(5);
  lights.reserve(5);
  materials.reserve(5);
  textures.reserve(3);

  // start with one light and material by default
  lights.push_back(LightInfo());
  materials.push_back(MaterialInfo());

  resetProjection(1.0f);
  resetView();
  // model is default identity
}

WorldState::~WorldState() {}

void WorldState::useProgram(unsigned int i)
{
  glUseProgram(shaders[i]);
  currentProgram = i;
}

void WorldState::setProgram(unsigned int i, GLuint program)
{
  while (shaders.size() <= i)
    shaders.push_back((GLuint)0);

  shaders[i] = program;
}

void WorldState::loadTransforms() { loadTransforms(model, view, projection); }

void WorldState::loadTransforms(glm::mat4 M)
{
  loadTransforms(M, view, projection);
}

void WorldState::loadTransforms(glm::mat4 M, glm::mat4 V, glm::mat4 P)
{
  // setup uniform variables
  glUniformMatrix4fv(
      glGetUniformLocation(shaders[currentProgram], "ModelMatrix"), 1, GL_FALSE,
      glm::value_ptr(M));
  glUniformMatrix4fv(
      glGetUniformLocation(shaders[currentProgram], "ModelViewMatrix"), 1,
      GL_FALSE, glm::value_ptr(V * M));
  glUniformMatrix4fv(
      glGetUniformLocation(shaders[currentProgram], "ProjectionMatrix"), 1,
      GL_FALSE, glm::value_ptr(P));
  glUniformMatrix4fv(glGetUniformLocation(shaders[currentProgram], "MVP"), 1,
                     GL_FALSE, glm::value_ptr(P * V * M));

  /* the Normal Matrix is the inverse transpose of the upper left 3x3 modelview
   * matrix
   * this is used to transform normals, as they do not transform the same way as
   * vertices. */
  glUniformMatrix3fv(
      glGetUniformLocation(shaders[currentProgram], "NormalMatrix"), 1, GL_TRUE,
      glm::value_ptr(glm::inverse(glm::mat3(V * M))));
}

void WorldState::loadObjectTransforms(glm::mat4 oM)
{
  loadTransforms(model * oM, view, projection);
}

void WorldState::loadLight(unsigned int i)
{
  // Load lights based on their index so as not to overwrite other lights.
  char lposition[20], la[20], ld[20], ls[20], intensity[20];
  sprintf(lposition, "Light%d.LPosition",
          i); // do not remove the \0 or you'll break shaders
  sprintf(la, "Light%d.La", i);
  sprintf(ld, "Light%d.Ld", i);
  sprintf(ls, "Light%d.Ls", i);
  sprintf(intensity, "Light%d.Intensity", i);

  glUniform4fv(glGetUniformLocation(shaders[currentProgram], lposition), 1,
               glm::value_ptr(lights[i].LPosition));
  glUniform3fv(glGetUniformLocation(shaders[currentProgram], la), 1,
               glm::value_ptr(lights[i].La));
  glUniform3fv(glGetUniformLocation(shaders[currentProgram], ld), 1,
               glm::value_ptr(lights[i].Ld));
  glUniform3fv(glGetUniformLocation(shaders[currentProgram], ls), 1,
               glm::value_ptr(lights[i].Ls));
  glUniform3fv(glGetUniformLocation(shaders[currentProgram], intensity), 1,
               glm::value_ptr(lights[i].Intensity));
}

void WorldState::loadLights()
{
  for (unsigned int i = 0; i < lights.size(); i++)
    loadLight(i);
}

void WorldState::loadMaterial(unsigned int i)
{
  // Load lights based on their index so as not to overwrite other lights.
  char ka[20], kd[20], ks[20], shininess[20];
  sprintf(ka, "Material%d.Ka", i);
  sprintf(kd, "Material%d.Kd", i);
  sprintf(ks, "Material%d.Ks", i);
  sprintf(shininess, "Material%d.Shininess", i);

  // printf("%f %f %f\n", materials[i].Ka);

  glUniform3fv(glGetUniformLocation(shaders[currentProgram], ka), 1,
               glm::value_ptr(materials[i].Ka));
  glUniform3fv(glGetUniformLocation(shaders[currentProgram], kd), 1,
               glm::value_ptr(materials[i].Kd));
  glUniform3fv(glGetUniformLocation(shaders[currentProgram], ks), 1,
               glm::value_ptr(materials[i].Ks));
  glUniform1f(glGetUniformLocation(shaders[currentProgram], shininess),
              materials[i].Shininess);
}

void WorldState::loadMaterials()
{
  // TODO: each load material currently overwrites itself, need to fix shaders
  // and loadMaterial(i)
  for (unsigned int i = 0; i < materials.size(); i++)
    loadMaterial(i);
}

void WorldState::loadColorMaterial(glm::vec4 color)
{
  glUniform3fv(glGetUniformLocation(shaders[currentProgram], "Material0.Ka"),
               1, glm::value_ptr(glm::vec3(color)));
  glUniform3fv(glGetUniformLocation(shaders[currentProgram], "Material0.Kd"),
               1, glm::value_ptr(glm::vec3(0.6, 0.6, 0.6)));
  glUniform3fv(glGetUniformLocation(shaders[currentProgram], "Material0.Ks"),
               1, glm::value_ptr(glm::vec3(0.6, 0.6, 0.6)));
  glUniform1f(
      glGetUniformLocation(shaders[currentProgram], "Material0.Shininess"),
      0.5f);
}

void WorldState::loadTexture(unsigned int i)
{
  char texsampler[20];
  sprintf(texsampler, "TexSampler%d", textures[i].m_tex_num);

  glUniform1i(glGetUniformLocation(shaders[currentProgram], texsampler),
              textures[i].m_tex_num);
}

void WorldState::loadTextures()
{
  for (unsigned int i = 0; i < textures.size(); i++)
    loadTexture(i);
}

void WorldState::resetProjection(float aspectRatio)
{
  projection = glm::perspective(glm::radians(50.0f), aspectRatio, 0.1f, 40.0f);
}

void WorldState::resetView()
{
  // Setup camera position/orientation.
  view = glm::lookAt(glm::vec3(0.0f, 2.5f, 10.0f), // eye
                     glm::vec3(0.0f, 0.0f, 0.0f),  // centre
                     glm::vec3(0.0f, 1.0f, 0.0f)   // up
                     );
}

void WorldState::resetModel() { model = glm::mat4(); }

void WorldState::updateView(float vTheta, float vPhi, float vDepth,
                            glm::vec3 vPan)
{
  // Setup camera position/orientation.
  glm::mat4 camera = glm::lookAt(vDepth * cameraOrigin,      // eye
                                 glm::vec3(0, 0, 0),         // centre
                                 glm::vec3(0.0f, 1.0f, 0.0f) // up
                                 );
  view = glm::rotate(camera, vPhi, glm::vec3(1, 0, 0));
  view = glm::rotate(view, vTheta, glm::vec3(0, 1, 0));
  view = glm::translate(view, glm::mat3(glm::inverse(view)) * vPan);
}
