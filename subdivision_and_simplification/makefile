CXX= g++
# debug version
CFLAGS= --std=c++11 -g -O0 -Wall -Wpedantic -DDEBUG -DLINUX
LFLAGS= 
# optimized version
CFLAGS_OPT= --std=c++11 -Wall -O3 -DNDEBUG -DLINUX
LFLAGS_OPT= 

L_DIR=
INC_DIR= -I/usr/include/eigen3/

SRC= $(wildcard *.cxx) 
SRC_E= $(wildcard exec_*.cxx)
SRC_O= $(filter-out $(SRC_E), $(SRC))
OBJ= $(SRC_O:%.cxx=%.o)
OBJ_OPT = $(SRC_O:%.cxx=%.o_opt)
EXEC= $(SRC_E:%.cxx=%)
EXEC_OPT= $(SRC_E:%.cxx=%_opt)

all: $(EXEC)

# DEBUG VERSION RULES
exec_simplify: exec_simplify.o $(OBJ)
	$(CXX) $(CFLAGS) $(INC_DIR) $(L_DIR) -o $@ $^ $(LFLAGS)

exec_subdiv: exec_subdiv.o $(OBJ)
	$(CXX) $(CFLAGS) $(INC_DIR) $(L_DIR) -o $@ $^ $(LFLAGS)

exec_cartel: exec_cartel.o $(OBJ)
	make -C Cartel lib METHOD=dbg
	$(CXX) $(CFLAGS) $(INC_DIR) $(L_DIR) -o $@ $^ $(LFLAGS) Cartel/Cartel_dbg.a $(shell pkg-config --static --libs glew glfw3 gl)

exec_glut: exec_glut.o $(OBJ)
	$(CXX) $(CFLAGS) $(INC_DIR) $(L_DIR) -o $@ $^ $(LFLAGS)  -lGL -lGLU -lglut

%.o: %.cxx $(wildcard *.hxx)
	$(CXX) $(CFLAGS) $(INC_DIR) -c -o $@ $< 

# OPTIMIZED VERSION RULES
exec_subdiv_opt: exec_subdiv.o_opt $(OBJ_OPT)
	$(CXX) $(CFLAGS_OPT) $(INC_DIR) $(L_DIR) -o $@ $^ $(LFLAGS_OPT)

exec_simplify_opt: exec_simplify.o_opt $(OBJ_OPT)
	$(CXX) $(CFLAGS_OPT) $(INC_DIR) $(L_DIR) -o $@ $^ $(LFLAGS_OPT)

exec_cartel_opt: exec_cartel.o_opt $(OBJ_OPT)
	make -C Cartel lib METHOD=opt
	$(CXX) $(CFLAGS) $(INC_DIR) $(L_DIR) -o $@ $^ $(LFLAGS) Cartel/Cartel_opt.a $(shell pkg-config --static --libs glew glfw3 gl)

exec_glut_opt: exec_glut.o_opt $(OBJ_OPT)
	$(CXX) $(CFLAGS) $(INC_DIR) $(L_DIR) -o $@ $^ $(LFLAGS)  -lGL -lGLU -lglut

%.o_opt: %.cxx $(wildcard *.hxx)
	$(CXX) $(CFLAGS_OPT) $(INC_DIR) -c -o $@ $< 

# MISC
TAGS: $(SRC)
	etags *.cxx *.hxx

clean:
	rm -rf *~ *.o *.o_opt *.obj *.vtk $(EXEC) $(EXEC_OPT) TAGS 

.PRECIOUS: $(OBJ) $(OBJ_OPT)

echo:
	@echo $(SRC)
	@echo $(SRC_E)
	@echo $(SRC_O)
	@echo $(EXEC)
