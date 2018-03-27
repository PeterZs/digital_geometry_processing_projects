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

#include "ShaderUtils.h"
#include "Utils.h"

#include <string>
#include <iostream>
#include <ctime>

void buildShader(GLenum shaderType, const char *path, GLuint &shader)
{
    GLint  length        = 0;
    GLint  result        = GL_FALSE;
    const GLchar *source = NULL;
    
    // create the gl shader symbol
    if (!shader)
    {
        shader = glCreateShader(shaderType);
        if (!shader)
        {
            printf("unable to create gl shader object\n");
            return;
        }
    }

    // load the source from file
    source = loadFileAsString(path, &length);
    if (!source)
    {
        glDeleteShader(shader);
        shader = 0;
        printf("unable to load shader from file\n");
        return;
    }

    // compile the shader
    glShaderSource(shader, 1, &source, &length);
    glCompileShader(shader);
    delete [] source;

    // if the compile failed, get the GL internal
    // debug log
    glGetShaderiv(shader, GL_COMPILE_STATUS, &result);
    if (result == GL_FALSE)
    {
        char *log;

        /* get the shader info log */
        glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &length);
        log = new char[length];
        glGetShaderInfoLog(shader, length, &result, log);

        /* print an error message and the info log */
        printf("Unable to compile the shader from %s: %s\n", path, log);
        delete [] log;

        glDeleteShader(shader);
        shader = 0;
        return;
    }

#ifdef DEBUG
    GLenum status = GL_NO_ERROR;
    if ((status = glGetError()) != GL_NO_ERROR)
        printf("Shader creation failed: %08x\n", status);
#endif
}

GLuint buildProgram(int num, GLuint *shaders)
{
    GLuint shaderProgram = glCreateProgram();
    if (!shaderProgram)
    {
        printf("Failed to create GL handles for pipeline, aborting\n");
        if (shaderProgram)
            glDeleteProgram(shaderProgram);
        return 0;
    }

    // Link the shader program
    for (int i = 0; i < num; i++)
    {
        if (shaders[i] == 0)
        {
            printf("Invalid shader, aborting\n");
            return 0;
        }
        glAttachShader(shaderProgram, shaders[i]);
    }
    glLinkProgram(shaderProgram);

    // ensure the linking worked
    GLint link_status;
    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &link_status);
    if (link_status == GL_FALSE)
    {
        char *log;
        GLint length;

        // get the program log info
        glGetProgramiv(shaderProgram, GL_INFO_LOG_LENGTH, &length);
        log = new char[length];
        glGetProgramInfoLog(shaderProgram, length, &link_status, log);

        /* print an error message and the info log */
        printf("Failing linking the shader program: %s\n", log);
        delete [] log;

        glDeleteProgram(shaderProgram);

        return 0;
    }

#ifdef DEBUG
    GLenum status = GL_NO_ERROR;
    if ((status = glGetError()) != GL_NO_ERROR)
        printf("Error creating shader program: %08x\n", status);
#endif

    return shaderProgram;
}
