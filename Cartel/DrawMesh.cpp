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

#include "DrawMesh.h"
#include "Utils.h"
#include "glm/glm.hpp"
#include <fstream>

DrawMesh::~DrawMesh() {}

void DrawMesh::init(int num_buffers)
{
  m_state->bindVAO();

  m_vbos.resize(num_buffers);

  // set the buffer IDs that the buffer objects will use
  m_ibo.setRender(0, *m_state);
  for (int i = 0; i < num_buffers; i++)
    m_vbos[i].setRender(i + 1, *m_state);
}

void DrawMesh::drawMesh(GLenum primitive)
{
  m_state->bindVAO();
  syncGPU(0, 0);
  glDrawElements(primitive, m_num_elem, GL_UNSIGNED_INT, 0);
}

void DrawMesh::addBuffer(int buff_num)
{
  m_vbos.push_back(VBuffer(buff_num, *m_state));
  m_vbos[buff_num].m_renderID = m_vbos[buff_num].m_renderID + 1;
}

void DrawMesh::loadVBuffer(unsigned int buffer_num, int data_size,
                           GLubyte *data, int data_offset, int num_attr,
                           attrib_info *attr_info)
{
  if (buffer_num >= m_vbos.size())
    addBuffer(buffer_num);

  m_vbos[buffer_num].loadBuffer(data_size, data, data_offset, num_attr,
                                attr_info);
}

void DrawMesh::loadIBuffer(int num_elem, int elem_size, int *data)
{
  // no need to bind vertex array buffer, since this is an index buffer
  m_ibo.loadBuffer(num_elem, elem_size, data);

  m_num_elem = num_elem;
}

// by having a range, we can choose to only update the buffers we have changed
// extent is the number of buffers after the base index to update
void DrawMesh::syncGPU(unsigned int base, size_t extent)
{
  if (base + extent > m_vbos.size() || extent == 0)
    extent = m_vbos.size() - base;

  // bind the vertex array object to store all the vertex settings
  m_state->bindVAO();
  for (unsigned int i = base; i < extent; i++)
    m_vbos[i].SyncBuffer();

  m_ibo.SyncBuffer();
}