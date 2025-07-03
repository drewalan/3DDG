#makefile
#makefile for 3D-DG
#Drew Murray
#11/29/16

###################USAGE######################
# to run, use command:
#	$make -f makefile


.PHONY: all clean
.SECONDARY: main-build

all: pre-build main-build post-build

pre-build:
	@echo

post-build:
	rm -f $(cppobj) $(hdf5obj)

main-build: 3D-DG

#HDF_INSTALL := /usr
#EXTLIB := -L$(HDF_INSTALL)/lib64
#INCLUDE := -I$(HDF_INSTALL)/include
#LIBSHDF := $(EXTLIB) $(HDF_INSTALL)/lib64/libhdf5.so

#C++ files
cppsrc = $(wildcard cppsrc/*.cpp)
cppsrc += $(wildcard cppsrc/*.h)
cppobj = $(wildcard cppsrc/*.o)
cppobj += $(wildcard *.o)
hdf5obj = $(wildcard cppsrc/*.h.gch)

#MPI_COMPILE_FLAGS = $(shell mpiCC.openmpi --showme:compile)
#MPI_LINK_FLAGS = $(shell mpiCC.openmpi --showme:link)

MPI_COMPILE_FLAGS = $(shell mpic++ --showme:compile)
MPI_LINK_FLAGS = $(shell mpic++ --showme:link)


#compiler settings/flags
CXXFLAGS = -std=c++0x -pg -isystem/usr/include/openmpi-x86_64 -L/usr/lib64/openmpi/lib -L/usr/lib64/openmpi/lib/openmpi   #debugging flags for valgrind
#CXXFLAGS = -std=c++0x -isystem/usr/include/openmpi-x86_64 -L/usr/lib64/openmpi/lib -L/usr/lib64/openmpi/lib/openmpi
#Warning flags, unused variable currently disabled b/c some functions are designed to be 'swappable' with others that need more arguments
CXXFLAGS += -Wall -Wextra -Wno-unused-but-set-variable -Wno-unused-parameter
CXX = h5c++
LDFLAGS = -lz -lm -lmpi #-lpng

3D-DG: $(cppsrc) makefile
	$(CXX) $(CXXFLAGS) $(MPI_COMPILE_FLAGS) $(cppsrc) $(LDFLAGS) $(MPI_LINK_FLAGS) -o $@
	@$(MAKE) --no-print-directory post-build
	

clean:
	rm -f $(cppobj) $(hdf5obj) 3D-DG
