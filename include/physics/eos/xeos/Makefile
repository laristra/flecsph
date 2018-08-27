##############################################################
#                                                            #
#                       GNU Makefile                         #
#                       OK 09.06.2017                        #
#                                                            #
##############################################################

#-- Path variables
ROOTDIR = .
OBJ_PATH = $(ROOTDIR)/obj
BIN_PATH = $(ROOTDIR)
SRC_PATH = $(ROOTDIR)/src
INC_PATH = $(ROOTDIR)/src
TEST_PATH = $(ROOTDIR)/tests
HEADERS = xeospoly.h xeoscoldn.h xeosbase.h physquan.h nvector.h
HEADERS+= filereader.h
INCHEADERS = pq_enum.inc quantities.inc
OBJECTS = $(addprefix $(OBJ_PATH)/, $(OBJS))
BINARIES = main.x
# # #
VPATH = $(CWD) $(SRC_PATH) $(OBJ_PATH) $(BIN_PATH) $(INC_PATH)

#-- Compiler, include directories and flags
INCS = -I$(INC_PATH)

##-- Intel compiler and flags
#CXX = icpc
#CXXFLAGS = -O2 -parallel   # optimizing
#CXXFLAGS = -g -O0          # debugging

#-- GNU compiler and flags: NOTE you need C++14-capable compiler!!!
CXX = g++-5
CXXFLAGS = -std=c++14 -O2

#-- Linking flags
LDFLAGS = -L$(OBJ_PATH) -L.

#-- Targets
.PHONY: all clean
all : obj main.x

obj : ; mkdir -p $(OBJ_PATH)

$(SRC_PATH)/physquan.h : $(INCHEADERS)

main.x : main.cc main.h $(HEADERS) $(INCHEADERS)
	$(CXX) -o $@ $(INCS) $(CXXFLAGS) $(LDFLAGS) main.cc

$(OBJ_PATH)/%.o : %.cc %.h
	$(CXX) -c -o $@ $(INCS) $(CXXFLAGS) $<

#-- Cleanup
clean:
	rm -rf $(BINARIES) $(OBJ_PATH)/*.o

