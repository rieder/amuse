ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Downloading the code
##### Remove if not needed #####
VERSION = 121257ede5486660d9eb5fbb23b2c377f02ae465

metisse.tar.gz:
	##### Modify URL as needed #####
	$(DOWNLOAD) https://github.com/TeamMETISSE/METISSE/archive/$(VERSION).tar.gz >$@

PATCHES := $(file < patches/series)
PATCHES := $(patsubst %,%,$(PATCHES))    # replace newlines with spaces

src/metisse-$(VERSION): metisse.tar.gz
	##### Modify as needed #####
	tar xf $<
	mv metisse-$(VERSION) src
	##### Apply patches here, if needed #####
	#patch <patch.patch
	#for p in $(PATCHES) ; do patch -p1 <patches/$$p ; done


# Building the code into a static library
DEPFLAGS += $(STOPCOND_CFLAGS) $(STOPCONDMPI_CFLAGS) $(AMUSE_MPI_CFLAGS)
DEPFLAGS += $(FORSOCKETS_CFLAGS) $(SIMPLE_HASH_CFLAGS) $(G6LIB_CFLAGS)
DEPFLAGS += $(OPENMP_FFLAGS)
FCFLAGS += $(DEPFLAGS) -ffpe-summary=all -fbacktrace -fcheck=all -g -Wall

LDLIBS += -lm $(STOPCOND_LIBS) $(STOPCONDMPI_LIBS) $(AMUSE_MPI_LIBS)
LDLIBS += $(FORSOCKETS_LIBS)

CODELIB := libmetisse.a

##### Modify build command as needed #####
.PHONY: src/METISSE/$(CODELIB)
src/METISSE/$(CODELIB): | src/METISSE/make
	$(MAKE) -C src/METISSE/make libmetisse


# Building the workers
metisse_worker.f90: interface.py
	amusifier --type=f90 interface.py MetisseInterface -o $@

storage.o: storage.f90
	$(MPIFC) -c -o $@ $(FCFLAGS) $<

metisse_worker.o: metisse_worker.f90 storage.o
	$(MPIFC) -c -o $@ $(FCFLAGS) $<

metisse_worker: interface.o metisse_worker.o storage.o src/METISSE/$(CODELIB)
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS)

%.o: %.f90
	$(MPIFC) -c -o $@ $(FCFLAGS) -Isrc/METISSE/make $<


# Which packages contain which workers?
amuse-metisse_contains: metisse_worker

# Building and installing packages
develop-%: %_contains
	support/shared/uninstall.sh $*
	python -m pip install -e packages/$*

install-%: %_contains
	support/shared/uninstall.sh $*
	python -m pip install packages/$*

package-%: %_contains
	python3 -m pip install -vv --no-cache-dir --no-deps --no-build-isolation --prefix ${PREFIX} packages/$*

test-%:
	cd packages/$* && pytest


# Cleaning up
.PHONY: clean
clean:
	$(MAKE) -C src clean
	rm -rf *.o *.mod *worker*

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache

