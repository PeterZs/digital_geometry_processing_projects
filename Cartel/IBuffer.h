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

#ifndef IBUFFER_H
#define IBUFFER_H

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "RenderState.h"

/**
 * The CPU side representation of a buffer used to store indices to index into
 * the other vertex attributes.
 * This has been made into a separate class for efficieny. It need not do as
 * much as a vertex buffer, and
 * is stored with different commands onto the GPU.
 */
class IBuffer
{
public:
  // takes in its index within the mesh that owns it
  // this is used to keep track of attrib array number
  IBuffer()
      : m_renderID(0), m_rState(NULL), m_local_data(NULL), m_elem_size(0),
        m_size(0)
  {
  }
  IBuffer(GLuint bufferID, RenderState &r_state)
      : m_renderID(bufferID), m_rState(&r_state), m_local_data(NULL),
        m_elem_size(0), m_size(0)
  {
  }
  ~IBuffer()
  {
    if (m_local_data)
      delete[] m_local_data;
  }

  void setRender(GLuint buffID, RenderState &r_state)
  {
    m_renderID = buffID;
    m_rState = &r_state;
  }
  void resizeBuffer(int size);
  void loadBuffer(int num_elem, int elem_size, int *data, int offset = 0);
  void SyncBuffer();
  void SyncBuffer(GLuint buffer);

  int m_size;            // the size in BYTEs of our local data store
  int m_elem_size;       // the size of an individual index
  GLubyte *m_local_data; // our local copy of the mesh data that we can modify
                         // and copy to the gpu

  GLuint m_renderID;
  RenderState *m_rState;
};

#endif // IBUFFER_H