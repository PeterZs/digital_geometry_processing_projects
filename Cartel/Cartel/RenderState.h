/* Copyright (c) Russell Gillette
 * December 2013
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

#pragma once
#ifndef RENDER_STATE_H
#define RENDER_STATE_H

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <vector>
using std::vector;


struct attrib_info
{
    unsigned int attrib_number;  // which attribute within the mesh this info is for
    unsigned int attrib_size;    // total size in bytes of one attribute
    unsigned int num_comp;       // the number of components each of these atributes holds
    unsigned int data_offset;    // the number of bytes from the start of the array to the first instance of this attribute
    unsigned int data_stride;    // the number of bytes between instances of this attribute
};

/**
 * A RenderState object stores a collection of GPU vertex buffer objects to be used to render a
 * single mesh at a time. It also stores all of the associated state setup for
 * those buffers, such as attribute stride and offset within a vertex
 * array object. 
 */
class RenderState
{
public:
    RenderState()
    {
        glGenVertexArrays(1, &m_vaoID);
    }
    ~RenderState()
    {
        glDeleteVertexArrays(1, &m_vaoID);

        if (!m_vboIDs.empty())
            glDeleteBuffers(m_vboIDs.size(), &m_vboIDs[0]);
    }
    void bindVAO()
    {
        glBindVertexArray(m_vaoID);
    }
    void bindVBO(unsigned int i, GLenum target = GL_ARRAY_BUFFER)
    {
        glBindBuffer(target, operator[](i));
    }

    void bindIBO(unsigned int i)
    { bindVBO(i, GL_ELEMENT_ARRAY_BUFFER); }


    void setNextBufferData(std::size_t num_bytes, unsigned char *data)
    {
        setBufferData(dirty++, num_bytes, data);
    }

    void setBufferData(unsigned int i, std::size_t num_bytes, unsigned char *data = NULL)
    {
        bindVBO(i);
        if (m_bytes[i] < num_bytes)
        {
            m_bytes[i] = num_bytes;
            glBufferData(GL_ARRAY_BUFFER, num_bytes, data, GL_STATIC_DRAW); // Static, because all changes happen in the edit mesh
        }
        else
            glBufferSubData(GL_ARRAY_BUFFER, 0, num_bytes, data);
    }
    void setAttributeData(attrib_info &info)
    { setAttributeData(info.attrib_number, info.num_comp, info.data_offset, info.data_stride); }

    void setAttributeData(unsigned int attrib_number, unsigned int num_comp, unsigned int data_offset, unsigned int data_stride)
    {
        glEnableVertexAttribArray(attrib_number);
        glVertexAttribPointer(attrib_number, num_comp, GL_FLOAT, GL_FALSE,
                              data_stride, (GLvoid *) data_offset);
    }

    GLuint operator[](unsigned int i)
    {
        if (i >= m_vboIDs.size())
        {
            m_vboIDs.resize(i+1, 0);
            m_bytes.resize(i+1, 0);
        }
        if (m_vboIDs[i] == 0)
            glGenBuffers(1, &m_vboIDs[i]);
        return m_vboIDs[i];
    }

    void makeClean()
    { dirty = 0; }

private:
    int dirty; // used for loading to copy into appropriate buffers

    GLuint m_vaoID;
    vector<GLuint> m_vboIDs;
    vector<GLuint> m_bytes; // used to avoid needless resizing of buffers
};

#endif // RENDER_STATE_H