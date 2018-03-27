/* Copyright (c) Russell Gillette
 * February 2014
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

#include "Texture.h"
#include <stdio.h>

int Texture::loadTexture(GLenum tex_num, texture_info &tex_info, void *data)
{
    return loadTexture(tex_num, GL_TEXTURE_2D, 0, tex_info, data);
}

int Texture::loadTexture(GLenum tex_num, GLuint target, texture_info &tex_info, void *data)
{
    return loadTexture(tex_num, target, 0, tex_info, data);
}

int Texture::loadTexture(GLenum tex_num, GLenum target, GLuint level, texture_info &tex_info, void *data)
{
    // TODO: check texture size before loading so as to
    //       not need to reallocate each time

    glActiveTexture(tex_num);
    m_tex_num = tex_num - GL_TEXTURE0;

    glEnable(m_target);
    glBindTexture(m_target, m_texID);
    m_tex_info = tex_info;

    // NOTE:
    // This sets all of the sample parameters that are built into textures.
    // Without doing this, tetures DO NOT WORK
    // There is an alternative, that binds this data outside of the texture for re-usability
    //    called a Sampler Object
    setDefaultSampleState();

    // populate our texture with the data from the Image to output our results to
    switch (target)
    {
    case GL_TEXTURE_1D_ARRAY:
    case GL_TEXTURE_2D:
    case GL_TEXTURE_CUBE_MAP_POSITIVE_X:
    case GL_TEXTURE_CUBE_MAP_NEGATIVE_X:
    case GL_TEXTURE_CUBE_MAP_POSITIVE_Y:
    case GL_TEXTURE_CUBE_MAP_NEGATIVE_Y:
    case GL_TEXTURE_CUBE_MAP_POSITIVE_Z:
    case GL_TEXTURE_CUBE_MAP_NEGATIVE_Z:
        {
            glTexImage2D(target, level, tex_info.int_format, tex_info.dimension[0], tex_info.dimension[1], tex_info.border, tex_info.ext_format, tex_info.type, data);
            break;
        }
    case GL_TEXTURE_CUBE_MAP:
        {
            printf("Failure, cannot map CUBE_MAP, must map individual faces\n");
            break;
        }
    // TODO:
    case GL_TEXTURE_1D:
    case GL_TEXTURE_2D_ARRAY:
    case GL_TEXTURE_3D:
    default:
        break;
    }

#ifdef DEBUG
    GLenum status = GL_NO_ERROR;
    if ((status = glGetError()) != GL_NO_ERROR)
        printf("Loading texture caused a gl error: %08x\n", status);
#endif

    return 1;
}

void Texture::setDefaultSampleState(GLenum wrap, GLenum filt)
{
    /***********************************************************
     * Set the Default Texture Sampler State
     * This is overridden if this texture is bound to a sampler
     ***********************************************************/
    // define textures to not repeat
    glTexParameteri( m_target, GL_TEXTURE_WRAP_S, wrap);
    glTexParameteri( m_target, GL_TEXTURE_WRAP_T, wrap);
    glTexParameteri( m_target, GL_TEXTURE_WRAP_R, wrap);
    // these are not technically needed, but will allow us to leverage hardware for
    // extra blurring by sampling between pixels (and using the linear interpolation)
    glTexParameteri( m_target, GL_TEXTURE_MIN_FILTER, filt);
    glTexParameteri( m_target, GL_TEXTURE_MAG_FILTER, filt);
}
