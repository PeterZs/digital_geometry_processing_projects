cmake_minimum_required(VERSION 2.6)

project( mmview )

## -------------------------------------------------------
## Debug
## -------------------------------------------------------

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif()

## -------------------------------------------------------
## GLUT and OpenGl
## -------------------------------------------------------

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

## -------------------------------------------------------
## Extra flags
## -------------------------------------------------------
include_directories(
  ${mmview_SOURCE_DIR}/mesh
  ${mmview_SOURCE_DIR}/mathvector
  ${OPENGL_INCLUDE_DIRS}
  ${GLUT_INCLUDE_DIRS}
)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${mmview_SOURCE_DIR})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_STANDARD 11)


## -------------------------------------------------------
## BUILD
## -------------------------------------------------------

SET(FILES_SOURCE
  ## glut
  glut/main.cpp
  glut/array3d.h
  glut/glmesh.h
  glut/glmesh_common.cpp
  glut/glmesh.cpp
  ## glmesh
  mathvector/mathvector.cpp
  mathvector/mathvector.h
  mathvector/matrix3.h
  ## mesh
  mesh/edge.cpp
  mesh/matrix3.h
  mesh/mesh.h
  mesh/mesh_io.h
  mesh/timeval.h
  mesh/triangle.h
  mesh/vertex.h
  mesh/edge.h
  mesh/mesh.cpp
  mesh/mesh_io.cpp
  mesh/misc.h
  mesh/triangle.cpp
  mesh/vertex.cpp
)

SET(TARGETLINKLIBS ${OPENGL_LIBRARIES}  ${GLUT_LIBRARIES} )
add_executable( mview ${FILES_SOURCE} )
target_link_libraries( mview ${TARGETLINKLIBS} )


## -------------------------------------------------------
## INSTALLING
## -------------------------------------------------------

##set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH})
##INCLUDE(InstallRequiredSystemLibraries)
## SET(CPACK_PACKAGE_NAME, "Cornucopia")
## SET(CPACK_PACKAGE_VENDOR "Ilya Baran")
## SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/gpl.txt")
## SET(CPACK_PACKAGE_VERSION_MAJOR "0")
## SET(CPACK_PACKAGE_VERSION_MINOR "1")
## SET(CPACK_PACKAGE_VERSION_PATCH "0")
## SET(CPACK_PACKAGE_EXECUTABLES "DemoUI" "Cornucopia DemoUI")
## INCLUDE( CPack )
