# --------------------------------------------------------------
#                              Petsc
# --------------------------------------------------------------

# Either set these here, or through environment variables
#PETSC_DIR=/home/hooshi/code/PETSc/petsc-3.6.3
#PETSC_ARCH=serial-opt
# Petsc 3.5 and newer
#include ${PETSC_DIR}/lib/petsc/conf/variables
# Petsc 3.4 and older
include ${PETSC_DIR}/conf/variables

# --------------------------------------------------------------
#                              Flags
# --------------------------------------------------------------

CC=g++
EXE=mview

ARMAFLAGS= -DARMA_DONT_USE_WRAPPER
#CPPFLAGS = $(INCPATH) -O2 $(ARMAFLAGS)
#MYCPPFLAGS = --std=c++11 -O0 -g3 -Wall $(ARMAFLAGS) $(INCPATH)  
MYCPPFLAGS = --std=c++11 ${CXX_FLAGS} $(ARMAFLAGS)  $(INCPATH) 
INCPATH = -I. -Imesh -Imathvector $(PETSC_CC_INCLUDES)
LIBS = -lGL -lGLU -lglut $(PETSC_LIB)

# --------------------------------------------------------------
#                              Rules
# --------------------------------------------------------------

GLUT =  glut/main.o \
	glut/glmesh.o \
	glut/glmesh_common.o 

MESH =	mesh/mesh.o \
	mesh/triangle.o \
	mesh/vertex.o \
	mesh/edge.o \
	mesh/mesh_io.o \
	mesh/arap.o

VEC =	mathvector/mathvector.o

all: GLUT

GLUT: $(GLUT) $(MESH) $(VEC)
	$(CC) $(MYCPPFLAGS) $?  $(LIBS) -o $(EXE)


%.o : %.cpp
	$(CC) $(MYCPPFLAGS) $(FLAGS) -c -o $@ $<

tidy:
	find . -name "*cpp" -exec touch -c {} \;

clean:
	rm -f $(EXE) $(QT3) $(QT3_MOC) $(QT4) $(QT4_MOC) $(GLUT) $(MESH) $(VEC)

echo:
	echo $(PETSC_CCPPFLAGS)
	echo ${CXX_FLAGS}
	echo ${CXXFLAGS}
	echo ${CCPPFLAGS}
