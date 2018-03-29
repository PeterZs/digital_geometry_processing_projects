/* Copyright (c) Russell Gillette
 * December 2013
 *
 * CPSC Computer Graphics 2014
 * University of British Columbia
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

#ifndef VBUFFER_H
#define VBUFFER_H

#include <cstring>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define MAX_NUM_ATTRIB 4

#include "RenderState.h"

/**
 * CPU side storage of the data to be stored in one VBO on the GPU. This object
 * can handle
 * many different vertex attribute layouts.
 */
class VBuffer
{
public:
  // we dont want the vbuffer owning the gl buffer id, beause then multiple
  // animations would ech require a different id (or deleting one would
  // inherently delete the gl id its working with)

  VBuffer()
      : m_renderID(0), m_rState(NULL), m_local_data(NULL), m_attr_info(NULL),
        m_num_attr(0), m_size(0)
  {
  }
  VBuffer(GLuint bufferID, RenderState &r_state)
      : m_renderID(bufferID), m_rState(&r_state), m_local_data(NULL),
        m_attr_info(NULL), m_num_attr(0), m_size(0)
  {
  }
  VBuffer(const VBuffer &old)
      : m_renderID(old.m_renderID), m_rState(old.m_rState),
        m_num_attr(old.m_num_attr), m_size(old.m_size)
  {
    m_local_data = (float *)new char[m_size];
    memcpy(m_local_data, old.m_local_data, m_size);

    m_attr_info = new attrib_info[m_num_attr];
    memcpy(m_attr_info, old.m_attr_info, sizeof(attrib_info) * m_num_attr);
  }

  ~VBuffer()
  {
    if (m_local_data)
      delete[] m_local_data;
    if (m_attr_info)
      delete[] m_attr_info;
  }

  void setRender(GLuint buffID, RenderState &r_state)
  {
    m_renderID = buffID;
    m_rState = &r_state;
  }

  void resizeBuffer(int size);
  // IMPORTANT: you must load the correct VAO prior to calling this function.
  // Since the
  // buffer object has no notion of VAO it is just setting its parameters on
  // whatever
  // VAO is currently bound
  void loadBuffer(int data_size, GLubyte *data, int data_offset, int num_attr,
                  attrib_info *attr_info);
  void SyncBuffer();
  void SyncBuffer(GLuint buffer);
  void SyncBuffer(GLuint buffer, GLuint *size);

  int m_size;     // the size in BYTEs of our local data store
  int m_num_attr; // the number of attributes per vertex
  attrib_info *m_attr_info;

  float *m_local_data; // our local copy of the mesh data that we can modify and
                       // copy to the gpu

  GLuint m_renderID;
  RenderState *m_rState;
};

#endif // VBUFFER_H
