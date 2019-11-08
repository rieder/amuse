# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=$(shell amusifier --get-amuse-dir)
-include $(AMUSE_DIR)/config.mk

MPICXX   ?= mpicxx

CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

PS_PATH = -Isrc/pentacle/fdps/src/
PPPT_PATH = -Isrc/pentacle/pppt/
use_x86 = yes
PHANTOM = -DUSE_PHANTOM

PS_PATH += -Isrc/pentacle/pppt/
PS_PATH += -Isrc/pentacle/hermite/
PS_PATH += -Isrc/pentacle/phantom/

CODE_PATH = -Isrc/pentacle/dens_system/

ORDER_MULTIPOLE = -DUSE_QUAD

ifeq ($(use_x86),yes)
CXX = time mpicxx
CXXFLAGS = -O3
#CXXFLAGS = -O2
#CXXFLAGS = -O1
CXXFLAGS += -Wall
#CXXFLAGS += -march=core-avx2
CXXFLAGS += -ffast-math -funroll-loops
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -std=c++14
CXXFLAGS += -Wall
#CXXFLAGS += -D_GLIBCXX_PARALLEL
#CXXFLAGS += -DINTRINSIC_X86
#CXXFLAGS += -DP3T_GPU_BARRIER_FOR_PROFILE
#CXXFLAGS += -DUSE_GNU_PARALLEL_SORT
endif

SRC = src/pentacle/dens_system/main.cpp

CODELIB = src/libpentacle.a

OBJS = interface.o

all: pentacle_worker 

clean:
	$(RM) -rf __pycache__
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ pentacle_worker worker_code.cc
	make -C src clean

distclean: clean
	make -C src distclean

$(CODELIB): $(SRC)
	$(MPICXX) $(PPPT_PATH) $(ORDER_MULTIPOLE) $(FORCE_PRECISION) $(PHANTOM) $(PS_PATH) $(CXXFLAGS) -c -o $@ $< $(CLIBS)

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py PentacleInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py PentacleInterface -o $@

pentacle_worker: worker_code.cc worker_code.h $(OBJS)
	$(MPICXX) $(PPPT_PATH) $(ORDER_MULTIPOLE) $(FORCE_PRECISION) $(PHANTOM) $(PS_PATH) $(CXXFLAGS) $< $(OBJS) -o $@

.cc.o: $<
	$(MPICXX) $(PPPT_PATH) $(CODE_PATH) $(ORDER_MULTIPOLE) $(FORCE_PRECISION) $(PHANTOM) $(PS_PATH) $(CXXFLAGS) -c -o $@ $< 
