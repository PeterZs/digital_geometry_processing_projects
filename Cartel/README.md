#Cartel
Mesh Viewing and Modification Tool

## Overview 
Cartel is a mesh visualization and modification tool designed for simplicity and clarity. It is built on
top of a basic OpenGL 3.3 environment using the GLFW windowing system, and enables complex
rendering effects through the easy use of shaders. Deployment and compilation are
simple due to the use of few, and common, dependencies.

## About
Cartel was developed in early 2014 by Russell Gillette and Darcy Harisson with the sole purpose of 
never having to use Meshlab ever again. It was developped from a rendering framework that Russell 
had been working on across many projects and a half-edge data structure that Darcy had whipped together 
to try and use Maya instead. Where meshlab offered an exceedingly complicated interface, exceedingly 
complicated dependancies, and minimal documentation, we have strived to offer a clean interface with 
only the basic modeling accoutrements. There is no QT, and model modification code is added directly to 
the model representation internally. Furthermore we enforce mesh manifoldness at all times. 
You will curse the asserts that this causes. But trust me, having lots of problems as soon as you do 
something stupid is better than having sparadic problems long after with no clue as to why. --Russell

This fork focuses on many usability improvements, both during setup, while coding, and 
while using the program:
* Improved Unix support with revamped Makefile. It *should* work on OSX just by installing few 
common dependencies. For more info, check the Makefile directly.
* Cleanup and code fixes
* Added ImGui as a user interface to improve usability such as loading models and changing view modes. 
The developer is also encouraged to add their own functionality in the menu bar, and can potentially 
implement controls such as sliders, inputs, etc to change parameters without modifying the source code.

--Giorgio

## Dependencies
The libraries included with the build are for VisualStudio 2012 or later, and the source code includes
some (though not many) C++ 11 operations, and thus will not work with earlier compilers. That
said, there is nothing OS specific within the code, and alternative libraries are easily acquired at the
below links.

### GLEW 
GLEW is the extension wrangler library for OpenGL that allows the use of newer OpenGL functionality
on computers that do not support it. You will never have to do anything using this library,
its just included, and thus worth mentioning. Its home webpage is here:
http://glew.sourceforge.net/

### GLM
GLM is the open gl math extensions and provides a lot of useful functionality when interfacing with
GLSL, the OpenGL shader language. This library provides types, and functions on those types
used SOLELY FOR RENDERING. This is important, because the glm types do not oer double
precision, and do not interface easily with Eigen. Use glm only for visual eects or when passing
data to the GPU in main.
http://glm.g-truc.net/0.9.5/index.html

### GLFW
GLFW is the windowing system used by Cartel. It is a common replacement for GLUT, and was
chosen due to its nicer rendering loop and newer interface. This windowing library is cross-platform,
allowing the deployment of this system across OSes, but does, at the time of writing this, make
setting up this project on a linux system with integrated graphics dificult.
http://www.glfw.org/

### Eigen
Eigen is the back-end linear-algebra library used by Cartel, and the largest source of complexity
within the project. This library makes heavy use of templates, but is well documented. Documentation
can found here: http://eigen.tuxfamily.org/index.php?title=Main_Page

### ImGui
ImGui is the user interface library used by Cartel. It should never create dependency issues as it's small 
and included in the project. For more information on the library check its page: https://github.com/ocornut/imgui
