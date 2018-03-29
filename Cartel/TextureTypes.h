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

#include "Texture.h"

#pragma once
#ifndef TEXTURE_TYPES_H
#define TEXTURE_TYPES_H

/* This class encapsulates data for a 2D texture. The goal with this class is
 * to construct a number of reasonable defaults and then allow the user to
 * modify only the data they need to via an easy API */
class Texture2D : public Texture
{
public:
  Texture2D() : Texture() {}
  Texture2D(GLuint tID) : Texture(tID, GL_TEXTURE_2D) {}

  /* load the texture to the texture slot tex_num, on target target
   * ex: texture slot: GL_TEXTURE0, target: GL_TEXTURE_2D */
  int loadTexture(GLenum tex_num, texture_info &tex_info, void *data)
  {
    int tx = Texture::loadTexture(tex_num, GL_TEXTURE_2D, 0, tex_info, data);
    setDefaultSampleState(GL_REPEAT, GL_NEAREST);
    return tx;
  }

  int loadTexture(GLenum tex_num, GLuint level, texture_info &tex_info,
                  void *data)
  {
    return Texture::loadTexture(tex_num, GL_TEXTURE_2D, level, tex_info, data);
  }

private:
};

/* This class encapsulates data for a 2D texture. The goal with this class is
 * to construct a number of reasonable defaults and then allow the user to
 * modify only the data they need to via an easy API */
class TextureCubeMap : public Texture
{
public:
  TextureCubeMap() : Texture(0, GL_TEXTURE_CUBE_MAP) {}
  TextureCubeMap(GLuint tID) : Texture(tID, GL_TEXTURE_CUBE_MAP) {}

  /* load the texture to the texture slot tex_num, on target target
   * ex: texture slot: GL_TEXTURE0, target: GL_TEXTURE_2D */
  int loadTexture0(GLenum tex_num, texture_info &tex_info, void *data)
  {
    return Texture::loadTexture(tex_num, GL_TEXTURE_CUBE_MAP_POSITIVE_X,
                                tex_info, data);
  }
  int loadTexture1(GLenum tex_num, texture_info &tex_info, void *data)
  {
    return Texture::loadTexture(tex_num, GL_TEXTURE_CUBE_MAP_NEGATIVE_X,
                                tex_info, data);
  }
  int loadTexture2(GLenum tex_num, texture_info &tex_info, void *data)
  {
    return Texture::loadTexture(tex_num, GL_TEXTURE_CUBE_MAP_POSITIVE_Y,
                                tex_info, data);
  }
  int loadTexture3(GLenum tex_num, texture_info &tex_info, void *data)
  {
    return Texture::loadTexture(tex_num, GL_TEXTURE_CUBE_MAP_NEGATIVE_Y,
                                tex_info, data);
  }
  int loadTexture4(GLenum tex_num, texture_info &tex_info, void *data)
  {
    return Texture::loadTexture(tex_num, GL_TEXTURE_CUBE_MAP_POSITIVE_Z,
                                tex_info, data);
  }
  int loadTexture5(GLenum tex_num, texture_info &tex_info, void *data)
  {
    return Texture::loadTexture(tex_num, GL_TEXTURE_CUBE_MAP_NEGATIVE_Z,
                                tex_info, data);
  }

  // NOTE:
  // ability to assign multiple levels to a cubemap has been omitted for now, I
  // don't think it will be widely used

private:
};

#endif // TEXTURE_TYPES_H
