# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../..
-include $(AMUSE_DIR)/config.mk

MPIFC ?= mpif90
FC      = $(MPIFC)

LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODELIBDIR = src/
CODELIB = src/libphantom.a

FCFLAGS += -L${CODELIBDIR} -lphantom

CODEDIR = src/phantom

CODE_GENERATOR = $(AMUSE_DIR)/build.py

SETUP = cluster
SYSTEM = gfortran

DOWNLOAD_FROM_WEB = $(PYTHON) ./download_http.py

OBJECTS = \
          ${CODEDIR}/build/boundary.o \
          ${CODEDIR}/build/centreofmass.o \
          ${CODEDIR}/build/checkoptions.o \
          ${CODEDIR}/build/checksetup.o \
          ${CODEDIR}/build/commons.o \
          ${CODEDIR}/build/config.o \
          ${CODEDIR}/build/coolfunc.o \
          ${CODEDIR}/build/cooling.o \
          ${CODEDIR}/build/damping.o \
          ${CODEDIR}/build/datafiles.o \
          ${CODEDIR}/build/dens.o \
          ${CODEDIR}/build/density_profiles.o \
          ${CODEDIR}/build/deriv.o \
          ${CODEDIR}/build/directsum.o \
          ${CODEDIR}/build/dtype_kdtree.o \
          ${CODEDIR}/build/dust.o \
          ${CODEDIR}/build/energies.o \
          ${CODEDIR}/build/eos.o \
          ${CODEDIR}/build/eos_helmholtz.o \
          ${CODEDIR}/build/eos_mesa.o \
          ${CODEDIR}/build/eos_mesa_microphysics.o \
          ${CODEDIR}/build/evolve.o \
          ${CODEDIR}/build/evwrite.o \
          ${CODEDIR}/build/extern_Bfield.o \
          ${CODEDIR}/build/extern_binary.o \
          ${CODEDIR}/build/extern_corotate.o \
          ${CODEDIR}/build/extern_gnewton.o \
          ${CODEDIR}/build/extern_gwinspiral.o \
          ${CODEDIR}/build/extern_lensethirring.o \
          ${CODEDIR}/build/extern_neutronstar.o \
          ${CODEDIR}/build/extern_prdrag.o \
          ${CODEDIR}/build/extern_spiral.o \
          ${CODEDIR}/build/extern_staticsine.o \
          ${CODEDIR}/build/externalforces.o \
          ${CODEDIR}/build/fastmath.o \
          ${CODEDIR}/build/force.o \
          ${CODEDIR}/build/fs_data.o \
          ${CODEDIR}/build/geometry.o \
          ${CODEDIR}/build/gitinfo.o \
          ${CODEDIR}/build/growth.o \
          ${CODEDIR}/build/h2chem.o \
          ${CODEDIR}/build/h2cooling.o \
          ${CODEDIR}/build/icosahedron.o \
          ${CODEDIR}/build/initial.o \
          ${CODEDIR}/build/io.o \
          ${CODEDIR}/build/kdtree.o \
          ${CODEDIR}/build/kernel_cubic.o \
          ${CODEDIR}/build/leastsquares.o \
          ${CODEDIR}/build/libphantom-evolve.o \
          ${CODEDIR}/build/libphantom-splash.o \
          ${CODEDIR}/build/libphantom.o \
          ${CODEDIR}/build/linklist_kdtree.o \
          ${CODEDIR}/build/lumin_nsdisc.o \
          ${CODEDIR}/build/memory.o \
          ${CODEDIR}/build/mf_write.o \
          ${CODEDIR}/build/mol_data.o \
          ${CODEDIR}/build/mpi_balance.o \
          ${CODEDIR}/build/mpi_dens.o \
          ${CODEDIR}/build/mpi_derivs.o \
          ${CODEDIR}/build/mpi_domain.o \
          ${CODEDIR}/build/mpi_force.o \
          ${CODEDIR}/build/mpi_utils.o \
          ${CODEDIR}/build/nicil.o \
          ${CODEDIR}/build/nicil_supplement.o \
          ${CODEDIR}/build/options.o \
          ${CODEDIR}/build/part.o \
          ${CODEDIR}/build/partinject.o \
          ${CODEDIR}/build/physcon.o \
          ${CODEDIR}/build/prompting.o \
          ${CODEDIR}/build/ptmass.o \
          ${CODEDIR}/build/quitdump.o \
          ${CODEDIR}/build/random.o \
          ${CODEDIR}/build/readwrite_dumps.o \
          ${CODEDIR}/build/readwrite_infile.o \
          ${CODEDIR}/build/set_binary.o \
          ${CODEDIR}/build/set_disc.o \
          ${CODEDIR}/build/set_dust.o \
          ${CODEDIR}/build/set_dust_options.o \
          ${CODEDIR}/build/set_flyby.o \
          ${CODEDIR}/build/set_slab.o \
          ${CODEDIR}/build/set_sphere.o \
          ${CODEDIR}/build/set_unifdis.o \
          ${CODEDIR}/build/set_vfield.o \
          ${CODEDIR}/build/setup_unifdis.o \
          ${CODEDIR}/build/solvelinearsystem.o \
          ${CODEDIR}/build/sort_particles.o \
          ${CODEDIR}/build/stack.o \
          ${CODEDIR}/build/step_leapfrog.o \
          ${CODEDIR}/build/step_supertimestep.o \
          ${CODEDIR}/build/stretchmap.o \
          ${CODEDIR}/build/test_cooling.o \
          ${CODEDIR}/build/test_corotate.o \
          ${CODEDIR}/build/test_derivs.o \
          ${CODEDIR}/build/test_dust.o \
          ${CODEDIR}/build/test_eos.o \
          ${CODEDIR}/build/test_externf.o \
          ${CODEDIR}/build/test_fastmath.o \
          ${CODEDIR}/build/test_geometry.o \
          ${CODEDIR}/build/test_gnewton.o \
          ${CODEDIR}/build/test_gravity.o \
          ${CODEDIR}/build/test_growth.o \
          ${CODEDIR}/build/test_indtstep.o \
          ${CODEDIR}/build/test_kdtree.o \
          ${CODEDIR}/build/test_kernel.o \
          ${CODEDIR}/build/test_link.o \
          ${CODEDIR}/build/test_luminosity.o \
          ${CODEDIR}/build/test_nonidealmhd.o \
          ${CODEDIR}/build/test_ptmass.o \
          ${CODEDIR}/build/test_rwdump.o \
          ${CODEDIR}/build/test_sedov.o \
          ${CODEDIR}/build/test_setdisc.o \
          ${CODEDIR}/build/test_step.o \
          ${CODEDIR}/build/testsuite.o \
          ${CODEDIR}/build/timestep.o \
          ${CODEDIR}/build/units.o \
          ${CODEDIR}/build/utils_allocate.o \
          ${CODEDIR}/build/utils_cpuinfo.o \
          ${CODEDIR}/build/utils_datafiles.o \
          ${CODEDIR}/build/utils_dumpfiles.o \
          ${CODEDIR}/build/utils_filenames.o \
          ${CODEDIR}/build/utils_indtimesteps.o \
          ${CODEDIR}/build/utils_infiles.o \
          ${CODEDIR}/build/utils_mathfunc.o \
          ${CODEDIR}/build/utils_omp.o \
          ${CODEDIR}/build/utils_sort.o \
          ${CODEDIR}/build/utils_sphNG.o \
          ${CODEDIR}/build/utils_spline.o \
          ${CODEDIR}/build/utils_summary.o \
          ${CODEDIR}/build/utils_supertimestep.o \
          ${CODEDIR}/build/utils_tables.o \
          ${CODEDIR}/build/utils_testsuite.o \
          ${CODEDIR}/build/utils_timing.o \
          ${CODEDIR}/build/utils_vectors.o \
          ${CODEDIR}/build/viscosity.o \
          ${CODEDIR}/build/writeheader.o
          # ${CODEDIR}/build/phantom.o \
          # ${CODEDIR}/build/getmathflags.o \



all: phantom_worker 

$(CODEDIR)/Makefile:
	make -C . download

download:
	$(RM) -Rf src
	mkdir src
	$(DOWNLOAD_FROM_WEB)

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ worker_code worker_code.f90
	$(RM) ${CODELIB}
	make -C ${CODEDIR} clean

$(CODELIB):$(CODEDIR)/Makefile
	make -C ${CODEDIR} libphantom SETUP=amuse SYSTEM=gfortran OPENMP=yes
	ar crs ${CODELIB} ${OBJECTS}

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py PhantomInterface -o $@

phantom_worker: worker_code.f90 $(CODELIB) $(OBJS)
	$(MPIFC) $(FCFLAGS) $(FS_FLAGS) $< $(OBJS) $(CODELIB) $(SC_FCLIBS) $(FS_LIBS) -fopenmp -o $@

%.o: %.f90
	$(FC) $(SC_FLAGS) $(FCFLAGS) -c -o $@ $<
