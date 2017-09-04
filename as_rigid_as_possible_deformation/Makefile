# --------------------------------------------------------------
#                              Flags
# --------------------------------------------------------------

CC=g++
EXE=mview

#CPPFLAGS = $(INCPATH) -O2 $(ARMAFLAGS)
#MYCPPFLAGS = --std=c++11 -O0 -g3 -Wall $(ARMAFLAGS) $(INCPATH)  
MYCPPFLAGS = --std=c++11 ${CXX_FLAGS}   $(INCPATH) 
INCPATH = -I. -Imesh -Imathvector 
LIBS = -lGL -lGLU -lglut 

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
	mesh/mesh_io.o

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
