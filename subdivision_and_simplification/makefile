CXX= g++
# debug version
CFLAGS= --std=c++11 -g3 -O0 -Wall -Wpedantic
LFLAGS= 
# optimized version
CFLAGS_OPT= --std=c++11 -g0 -O2 -DNDEBUG
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
exec_%: exec_%.o $(OBJ)
	$(CXX) $(CFLAGS) $(INC_DIR) $(L_DIR) -o $@ $^ $(LFLAGS)

%.o: %.cxx $(wildcard *.hxx)
	$(CXX) $(CFLAGS) $(INC_DIR) -c -o $@ $< 

# OPTIMIZED VERSION RULES
exec_%_opt: exec_%.o_opt $(OBJ_OPT)
	$(CXX) $(CFLAGS_OPT) $(INC_DIR) $(L_DIR) -o $@ $^ $(LFLAGS_OPT)

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
