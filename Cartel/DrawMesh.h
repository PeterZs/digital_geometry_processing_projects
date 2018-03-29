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
#pragma once
#ifndef DRAW_MESH_H
#define DRAW_MESH_H

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "IBuffer.h"
#include "RenderState.h"
#include "VBuffer.h"
#include <glm/glm.hpp>
#include <string>
#include <vector>

using namespace std;

/**
 * Storage for all of the vertex attributes needed for rending a single mesh,
 * as well as any state data specific to that mesh. (such as object model
 * transforms)
 */
class DrawMesh
{
public:
  glm::mat4 m_MV; // object's model transformations

  DrawMesh(RenderState &r_state) : m_num_elem(0), m_state(&r_state), m_vbos(0)
  {
  }
  ~DrawMesh();
  // takes the number of vbos and then a state which has enough space allocated
  // for those vbos
  // AND an ibo
  void init(int num_buffers);
  void drawMesh(GLenum primitive = GL_TRIANGLES);

  void addBuffer(int buff_num);

  // these functions copy the the data to our internally managed memory, to
  // modify
  // existing memory, use "GetData()
  void loadVBuffer(
      unsigned int buffer_num, int size, GLubyte *data, int data_offset,
      int num_attr,
      attrib_info *attr_info); // load data into the specified vertex buffer
  void loadIBuffer(int num_elem, int elem_size,
                   int *data); // load data into the index buffer

  // return a pointer to our internal copy of buffer i
  const float *getData(int i) const { return m_vbos[i].m_local_data; }
  float *getData(int i) { return m_vbos[i].m_local_data; }

  int getNumVBO() { return m_vbos.size(); }
  void
  syncGPU(unsigned int base,
          size_t extent = 0); // copy the changes to our local data to the GPU

private:
  void resizeBuffer(int size); // resize the internal store for the buffer

  GLsizei m_num_elem;
  IBuffer m_ibo;

  RenderState *m_state;
  vector<VBuffer> m_vbos;
};

#endif // DRAW_MESH