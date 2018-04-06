# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include $(AMUSE_DIR)/config.mk

MPICXX = mpicxx
MPICC = mpicc
CC = $(MPICC)
CXX = $(MPICXX)

CFLAGS   += -g #-Wall -g
LDFLAGS  += -lm $(MUSE_LD_FLAGS)


OBJS = interface.o

CODEDIR = src/gasoline
MDLDIR = src/mdl/null

CODE_GENERATOR = $(AMUSE_DIR)/build.py

DOWNLOAD_FROM_WEB = $(PYTHON) ./download_http.py


#	  /$$$$$$                                /$$ /$$                    
#	 /$$__  $$                              | $$|__/                    
#	| $$  \__/  /$$$$$$   /$$$$$$$  /$$$$$$ | $$ /$$ /$$$$$$$   /$$$$$$ 
#	| $$ /$$$$ |____  $$ /$$_____/ /$$__  $$| $$| $$| $$__  $$ /$$__  $$
#	| $$|_  $$  /$$$$$$$|  $$$$$$ | $$  \ $$| $$| $$| $$  \ $$| $$$$$$$$
#	| $$  \ $$ /$$__  $$ \____  $$| $$  | $$| $$| $$| $$  | $$| $$_____/
#	|  $$$$$$/|  $$$$$$$ /$$$$$$$/|  $$$$$$/| $$| $$| $$  | $$|  $$$$$$$
#	 \______/  \_______/|_______/  \______/ |__/|__/|__/  |__/ \_______/
                                                                                     
# Makefile for Gasoline (specify below)
# Then type make <platform>
# <platform> should be one of:
# null, pvm, pthread, ampi, lam_mpi, qmpi, mpi,
# charm
#

# Which C compiler should I use?
#CC = mpicc
#CC = gcc

#	 /$$      /$$                 /$$          
#	| $$$    /$$$                | $$          
#	| $$$$  /$$$$  /$$$$$$   /$$$$$$$  /$$$$$$ 
#	| $$ $$/$$ $$ /$$__  $$ /$$__  $$ /$$__  $$
#	| $$  $$$| $$| $$  \ $$| $$  | $$| $$$$$$$$
#	| $$\  $ | $$| $$  | $$| $$  | $$| $$_____/
#	| $$ \/  | $$|  $$$$$$/|  $$$$$$$|  $$$$$$$
#	|__/     |__/ \______/  \_______/ \_______/
	                                           
# What kind of problem are we running?  This will enable a default set of code components
# appropriate for different kinds of problems.  Select only one.
MODE_DEF = -DGALAXY            #Default mode: includes SPH, star formation, and feedback.
#MODE_DEF = -DCOLLISIONAL		#Collisional mode for simulating asteroids & planets
#MODE_DEF = -DNBODY 			#N-body only mode: gravity with nothing else.
#MODE_DEF = -DGLASS_IC 			#Generate a gravity glass (useful for making ICS)

#	 /$$   /$$                                         /$$
#	| $$  /$$/                                        | $$
#	| $$ /$$/   /$$$$$$   /$$$$$$  /$$$$$$$   /$$$$$$ | $$
#	| $$$$$/   /$$__  $$ /$$__  $$| $$__  $$ /$$__  $$| $$
#	| $$  $$  | $$$$$$$$| $$  \__/| $$  \ $$| $$$$$$$$| $$
#	| $$\  $$ | $$_____/| $$      | $$  | $$| $$_____/| $$
#	| $$ \  $$|  $$$$$$$| $$      | $$  | $$|  $$$$$$$| $$
#	|__/  \__/ \_______/|__/      |__/  |__/ \_______/|__/

# Which SPH Kernel should we use?
KERNEL_DEF = -DWENDLANDC4 #Use the Wendland C_4 Kernel (See Dehnen & Aly 2012)
#KERNEL_DEF = -DWENDLAND              #Use the Wendland C_2 Kernel (See Dehnen & Aly 2012)
#KERNEL_DEF = -DHSHRINK               #M4 kernel times (pi/6)^(1/3) for a tighter kernel
#KERNEL_DEF = -DM43D                  #Use a fancier 3D derived M4 Kernel
#KERNEL_DEF = -DPEAKEDKERNEL          #Use a modified Peaked M4 kernel as per Thomas and Couchman 92 (Bibcode:1992MNRAS.257...11T)
#KERNEL_DEF = -DQUINTIC               #Use a quintic spline kernel from Dehnen & Aly 2012 (Bibcode:2012MNRAS.425.1068D)


#	  /$$$$$$                      /$$ /$$                    
#	 /$$__  $$                    | $$|__/                    
#	| $$  \__/  /$$$$$$   /$$$$$$ | $$ /$$ /$$$$$$$   /$$$$$$ 
#	| $$       /$$__  $$ /$$__  $$| $$| $$| $$__  $$ /$$__  $$
#	| $$      | $$  \ $$| $$  \ $$| $$| $$| $$  \ $$| $$  \ $$
#	| $$    $$| $$  | $$| $$  | $$| $$| $$| $$  | $$| $$  | $$
#	|  $$$$$$/|  $$$$$$/|  $$$$$$/| $$| $$| $$  | $$|  $$$$$$$
#	 \______/  \______/  \______/ |__/|__/|__/  |__/ \____  $$
#	                                                 /$$  \ $$
#	                                                |  $$$$$$/
#	                                                 \______/ 
#
# The following pairs of lines define which cooling method to use.  Uncomment
# only one pair

COOLING_OBJ = $(CODEDIR)/cooling_metal.o
COOLING_DEF = -DCOOLING_METAL

#COOLING_OBJ = cooling_cosmo.o
#COOLING_DEF = -DCOOLING_COSMO

#COOLING_OBJ = cooling_poly.o
#COOLING_DEF = -DCOOLING_POLY

#COOLING_OBJ = cooling_grackle.o
#COOLING_DEF = -DCOOLING_GRACKLE
#COOLING_LIB = libgrackle.so -lhdf5

#COOLING_OBJ = cooling_metal_H2.o
#COOLING_DEF = -DCOOLING_MOLECULARH

#COOLING_OBJ = cooling_planet.o
#COOLING_DEF = -DCOOLING_PLANET

#	 /$$$$$$$$             /$$                       
#	| $$_____/            | $$                       
#	| $$       /$$   /$$ /$$$$$$    /$$$$$$  /$$$$$$ 
#	| $$$$$   |  $$ /$$/|_  $$_/   /$$__  $$|____  $$
#	| $$__/    \  $$$$/   | $$    | $$  \__/ /$$$$$$$
#	| $$        >$$  $$   | $$ /$$| $$      /$$__  $$
#	| $$$$$$$$ /$$/\  $$  |  $$$$/| $$     |  $$$$$$$
#	|________/|__/  \__/   \___/  |__/      \_______/
#	                                                 
#	 /$$$$$$$$                    /$$                                            
#	| $$_____/                   | $$                                            
#	| $$     /$$$$$$   /$$$$$$  /$$$$$$   /$$   /$$  /$$$$$$   /$$$$$$   /$$$$$$$
#	| $$$$$ /$$__  $$ |____  $$|_  $$_/  | $$  | $$ /$$__  $$ /$$__  $$ /$$_____/
#	| $$__/| $$$$$$$$  /$$$$$$$  | $$    | $$  | $$| $$  \__/| $$$$$$$$|  $$$$$$ 
#	| $$   | $$_____/ /$$__  $$  | $$ /$$| $$  | $$| $$      | $$_____/ \____  $$
#	| $$   |  $$$$$$$|  $$$$$$$  |  $$$$/|  $$$$$$/| $$      |  $$$$$$$ /$$$$$$$/
#	|__/    \_______/ \_______/   \___/   \______/ |__/       \_______/|_______/ 
                                                                                                                              
EXTRA_DEF =                                                                  
## GRAVITY AND HYDRODYNAMICS (REQUIRE -DGASOLINE)
#EXTRA_DEF += -DCULLENDEHNEN          #Use the time-variable Cullen & Dehnen 2010 viscosity limit switch
#EXTRA_DEF += -DVARALPHA              #Use a Morris & Monaghan-style time variable viscosity
#EXTRA_DEF += -DINFLOWOUTFLOW         #Enable inflow/outflow boundary conditions
#EXTRA_DEF += -DOUTURBDRIVER          #Drive turbulence as in Price & Federrath 2010 with OU variables
#EXTRA_DEF += -DROT_FRAME             #Use a rotating reference frame with rotation rate given by dOmega
#EXTRA_DEF += -DSIMPLE_GAS_DRAG       #Apply a drag to the gas velocities using either the Epstein regime, or stopping time
#EXTRA_DEF += -DPEXT                  #Allow non-zero external pressure 
#EXTRA_DEF += -DJEANSSMOOTH 		  #Calculate presure floor to keep jeans length larger than SPH h only (default is larger of SPH h or eps)
#EXTRA_DEF += -DJEANSSOFTONLY         #Calculate presure floor to keep jeans length larger than eps only
#EXTRA_DEF += -DSINKINGAVERAGE        #Use a smoothed relative velocity rather than per-particle for SMBH accretion rates

## STAR FORMATION AND FEEDBACK (REQUIRE -DSTARFORM)
#EXTRA_DEF += -DBLASTWAVE			  #Use the Stinson+ 2006 Blastwave feedback 
#EXTRA_DEF += -DNONTHERMAL			  #Use an Agertz+ 2013/Teyssier+ 2013 non-thermal pressure feedback
#EXTRA_DEF += -DKROUPA01              #Use the newer Kroupa 01 IMF (DOI:10.1046/j.1365-8711.2001.04022.x)
#EXTRA_DEF += -DKROUPA                #Use the Kroupa IMF (from Kroupa et al 1993 Bibcode:1993MNRAS.262..545K)
#EXTRA_DEF += -DMILLERSCALO			  #Use the old-school Miller & Scalo 1979 IMF (DOI:10.1086/190629)
#EXTRA_DEF += -DFEEDBACKDIFFLIMIT     #Disable thermal diffusion for particles that have their cooling shut off (requires -DBLASTWAVE)
#EXTRA_DEF += -DTOPHATFEEDBACK        #Use a tophat kernel to smooth feedback over (default with superbubble feedback)
#EXTRA_DEF += -DVOLUMEFEEDBACK        #Smooth FB energy with a volume-weighted kernel
#EXTRA_DEF += -DMASSFEEDBACK          #Smooth FB energy over mass-weighted kernel (default with blastwave feedback)

## COLLISIONS AND PLANET FORMATION (REQUIRE -DCOLLISIONS)
#EXTRA_DEF += -DAGGS                  #Include support for aggregates (like asteroids and rubble piles) See aggs.c
#EXTRA_DEF += -DGROWMASS              #Allow collisional particles to grow
#EXTRA_DEF += -DN_DIM                 #Number of dimensions to use for solar system calculations (Defaults to 3)
#EXTRA_DEF += -DSURFACEAREA           #Calculate and output the surface area of all the particles.
#EXTRA_DEF += -DNORMAL                #Calculate and output surface normals for particles (Needs -DSURFACEAREA)
#EXTRA_DEF += -DRUBBLE_ZML            #Use "Rubble Pile" collision model 
#EXTRA_DEF += -DSAND_PILE             #Use the "Sand Pile" collision model 
#EXTRA_DEF += -DSLIDING_PATCH         #Alternative collision scheme (I don't really know what this does) 
#EXTRA_DEF += -DTUMBLER               #Generate a hard-walled cylinder boundary (for use with collisions and sand piles)

## OUTPUT VISUALIZATION
#EXTRA_DEF += -DGSS_DUMPFRAME         #Dumpframe that allows coloring by particle property (Greg's version)
#EXTRA_DEF += -DVOXEL                 #Dump voxels for volume imaging.

## INITIAL CONDITION BUILDING
#EXTRA_DEF += -DGLASSZ                #Build a glass only in the z-direction, useful for building disc initial conditions

## DEBUG
#EXTRA_DEF += -DSETTRAPFPE            #Enable Floating point exceptions

#  /$$       /$$ /$$                                    /$$                    
#  | $$      |__/| $$                                   |__/                    
#  | $$       /$$| $$$$$$$   /$$$$$$  /$$$$$$   /$$$$$$  /$$  /$$$$$$   /$$$$$$$
#  | $$      | $$| $$__  $$ /$$__  $$|____  $$ /$$__  $$| $$ /$$__  $$ /$$_____/
#  | $$      | $$| $$  \ $$| $$  \__/ /$$$$$$$| $$  \__/| $$| $$$$$$$$|  $$$$$$ 
#  | $$      | $$| $$  | $$| $$      /$$__  $$| $$      | $$| $$_____/ \____  $$
#  | $$$$$$$$| $$| $$$$$$$/| $$     |  $$$$$$$| $$      | $$|  $$$$$$$ /$$$$$$$/
#  |________/|__/|_______/ |__/      \_______/|__/      |__/ \_______/|_______/ 
#                                                                               

# If you do dump frame stuff with png you need these libs ( -DUSE_PNG )
#PNG_INCL = -I/usr/include
#PNG_LIB = -L/usr/lib -lpng -lz
#PNG_LIB = -lpng -lz
#PNG_OBJ = writepng.o
#PNG_DEF = $(PNG_INCL) -DUSE_PNG
PNG_INCL =
PNG_LIB =
PNG_OBJ =
PNG_DEF =

# if you use OUTURBDRIVER or TURBULENT you need libgsl
#TURB_OBJ = outurb.o
#GSL_LIB = -lgsl -lgslcblas
TURB_OBJ =
GSL_LIB =

BASE_LD_FLAGS = $(PNG_LIB) $(GSL_LIB)
BASE_LD_FLAGS = $(PNG_LIB) $(GSL_LIB) $(COOLING_LIB)

# Remove the object for NBODY runs
ifeq ($(strip $(MODE_DEF)),-DNBODY)
COOLING_OBJ =
endif

BASE_DEF = $(PNG_DEF) $(COOLING_DEF) $(MODE_DEF) $(KERNEL_DEF) $(EXTRA_DEF)

# The filename of the compiled library
CODELIB = src/libgasoline.a


#       Basic MPI defines
#
#MPI_MDL			= $(MDLDIR)/mpi
MDL_CFLAGS  = -I$(MDLDIR)
CFLAGS		+= -O3 $(BASE_DEF)
LD_FLAGS	+= $(BASE_LD_FLAGS)
XOBJ		+= $(CODEDIR)/v_sqrt1.o
LIBMDL		= mdl.o -lm
MPI_MDL_CFLAGS	= -O3


OBJ	= $(CODEDIR)/main.o $(CODEDIR)/master.o $(CODEDIR)/param.o \
      $(CODEDIR)/outtype.o $(CODEDIR)/pkd.o $(CODEDIR)/pst.o $(CODEDIR)/grav.o \
	  $(CODEDIR)/ewald.o $(CODEDIR)/walk.o $(CODEDIR)/eccanom.o \
      $(CODEDIR)/fdl.o $(CODEDIR)/htable.o $(CODEDIR)/smooth.o \
	  $(CODEDIR)/smoothfcn.o $(CODEDIR)/collision.o $(CODEDIR)/qqsmooth.o \
	  $(COOLING_OBJ) \
      $(CODEDIR)/cosmo.o $(CODEDIR)/romberg.o \
	  $(CODEDIR)/starform.o $(CODEDIR)/feedback.o $(CODEDIR)/millerscalo.o \
	  $(CODEDIR)/supernova.o $(CODEDIR)/supernovaia.o \
	  $(CODEDIR)/startime.o $(CODEDIR)/stiff.o $(CODEDIR)/runge.o \
	  $(CODEDIR)/dumpframe.o $(CODEDIR)/dffuncs.o $(CODEDIR)/dumpvoxel.o \
	  $(CODEDIR)/rotbar.o $(CODEDIR)/special.o $(CODEDIR)/ssio.o \
	  $(PNG_OBJ) $(TURB_OBJ) \
	  $(CODEDIR)/treezip.o $(CODEDIR)/log.o $(CODEDIR)/linalg.o $(CODEDIR)/aggs.o
#OBJ	= main.o master.o param.o outtype.o pkd.o pst.o grav.o \
#	  ewald.o walk.o eccanom.o fdl.o htable.o smooth.o \
#	  smoothfcn.o collision.o qqsmooth.o $(COOLING_OBJ) cosmo.o romberg.o \
#	  starform.o feedback.o millerscalo.o supernova.o supernovaia.o \
#	  startime.o stiff.o runge.o dumpframe.o dffuncs.o dumpvoxel.o \
#	  rotbar.o special.o ssio.o $(PNG_OBJ) $(TURB_OBJ) \
#	  treezip.o log.o linalg.o aggs.o

#EXTRA_OBJ = erf.o v_sqrt1.o 

all: gasoline_worker

#gasoline_worker: lib
gasoline_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CFLAGS) $(SC_FLAGS) $(LDFLAGS) -I$(CODEDIR) $< $(OBJS) $(CODELIB) $(SC_CLIBS) -o $@

$(CODELIB): mdl.o $(OBJ) $(XOBJ)
	ar ruv $(CODELIB) mdl.o $(OBJ) $(XOBJ)
	ranlib $(CODELIB)

mdl.o:# mdl.c mdl.h
	$(MPICC) $(CFLAGS) $(MDL_CFLAGS) -c $(MDLDIR)/mdl.c -o $@

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py GasolineInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=h interface.py GasolineInterface -o $@

clean:
	$(RM) -r -f *.so *.dSYM *.o *.pyc worker_code.cc interface.h
	$(RM) -f src/libgasoline.a
	$(RM) -f *~ gasoline_worker
	$(RM) -f $(OBJ) $(XOBJ) mdl.o

depend:
	makedepend -Y -- $(BASE_DEF) -- *.c

$(XDIR):
	-mkdir $(BDIR)
	-mkdir $(XDIR)

$(OBJ) $(XOBJ): Makefile

# DO NOT DELETE
cosmo.o: runge.h cosmo.h
ewald.o: ewald.h pkd.h floattype.h cooling.h meval.h qeval.h
fdl.o: htable.h fdl.h
grav.o: pkd.h floattype.h cooling.h grav.h meval.h qeval.h
htable.o: htable.h
integration.o: linalg.h floattype.h integration.h
aggs.o: aggs.c aggs.h
linalg.o: linalg.h floattype.h
main.o: master.h param.h pst.h pkd.h floattype.h cooling.h smoothfcn.h
main.o: parameters.h cosmo.h outtype.h
master.o: master.h param.h pst.h pkd.h floattype.h cooling.h smoothfcn.h
master.o: parameters.h cosmo.h tipsydefs.h opentype.h fdl.h htable.h
master.o: outtype.h
millerscalo.o: millerscalo.h
outtype.o: pkd.h floattype.h cooling.h outtype.h
outurb.o: pkd.h outurb.h
param.o: param.h
pkd.o: pkd.h floattype.h cooling.h ewald.h grav.h walk.h opentype.h
pkd.o: tipsydefs.h dumpframe.h
pst.o: pst.h pkd.h floattype.h cooling.h smoothfcn.h starform.h feedback.h
pst.o: outtype.h smooth.h dumpframe.h
romberg.o: floattype.h
smooth.o: smooth.h pkd.h floattype.h cooling.h smoothfcn.h
smoothfcn.o: smoothfcn.h pkd.h floattype.h cooling.h
supernova.o: pkd.h floattype.h cooling.h feedback.h supernova.h startime.h
supernova.o: millerscalo.h supernovaia.h
supernovaia.o: supernova.h startime.h millerscalo.h feedback.h pkd.h
supernovaia.o: floattype.h cooling.h supernovaia.h
dumpframe.o: dumpframe.c dumpframe.h
dumpvoxel.o: dumpvoxel.c dumpframe.h dumpvoxel.h
dffuncs.o: dffuncs.c dumpframe.h
treezip.o: treezip.c treezip.h treezipkey.h treeziptypes.h
writepng.o: writepng.c writepng.h
walk.o: walk.h pkd.h floattype.h cooling.h
rotbar.o: master.h rotbar.h
cooling_metal.o: cooling_metal.h cooling_metal.c
cooling_grackle.o: cooling_grackle.h cooling_grackle.c
log.o: log.h

interface.o: interface.c
	$(MPICC) $(CFLAGS) $(MDL_CFLAGS) -I$(CODEDIR) -c -o $@ $< 
