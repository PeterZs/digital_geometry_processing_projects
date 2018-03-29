/* Copyright (c) Russell Gillette
 * February 2014
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
#ifndef TEXTURE_H
#define TEXTURE_H

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

struct texture_info
{
  // the internal and externally exposed formats of this texture
  int int_format;
  int ext_format;

  int dimension[3]; // { width, height, depth } 0 if unused
  int border; // this is a rather unused feature, and for now is not exposed on
              // this API
  int type; // this specifies the data type for each component (eg. half float)

  texture_info()
      : int_format(GL_RGB), ext_format(GL_RGB), border(0),
        type(GL_UNSIGNED_BYTE)
  {
  }
  texture_info(int width, int height)
      : int_format(GL_RGB), ext_format(GL_RGB), border(0),
        type(GL_UNSIGNED_BYTE)
  {
    dimension[0] = width;
    dimension[1] = height;
  }
};

enum TEXTURE_FLAGS
{
  TEXFLAGS_NONE = 0x0,
  TEXFLAGS_LOCAL = 0x1, // store a local copy of the texture to be uploaded at
                        // user discretion
  TEXFLAGS_ALL = 0xffffffff,
};

/* The texture class is responsible for holding all fo the data for a texture,
 * where OpenGL's concept of a texture consists of all mip levels for the
 * texture,
 * and all components if it is a texture array or cube_map */
class Texture
{
public:
  GLenum m_tex_num;

  Texture() : m_texID(0), m_target(GL_TEXTURE_2D), m_tex_num(GL_TEXTURE0)
  {
    glGenTextures(1, &m_texID);
  }
  Texture(GLuint tID, GLenum target)
      : m_target(target), m_texID(tID), m_tex_num(GL_TEXTURE0)
  {
  }
  virtual ~Texture() { glDeleteTextures(1, &m_texID); }

  int setSampler();

  /* load the texture to the texture slot tex_num, on target target
   * ex: texture slot: GL_TEXTURE0, target: GL_TEXTURE_2D, if not
   * specified, the target defaults to GL_TEXTURE_2D, and level to 0 */
  int loadTexture(GLenum tex_num, texture_info &tex_info, void *data);
  int loadTexture(GLenum tex_num, GLuint target, texture_info &tex_info,
                  void *data);
  int loadTexture(GLenum tex_num, GLenum target, GLuint level,
                  texture_info &tex_info, void *data);

protected:
  GLenum m_target;
  GLuint m_texID;
  texture_info m_tex_info;

  void setDefaultSampleState(GLenum wrap = GL_CLAMP_TO_EDGE,
                             GLenum filt = GL_LINEAR);
};

#endif // TEXTURE_H