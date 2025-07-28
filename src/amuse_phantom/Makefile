ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Downloading the code
VERSION = 4c7499e2641b6d59b790dc71c30e819515bd4bfe

phantom.tar.gz:
	$(DOWNLOAD) https://github.com/rieder/phantom/archive/$(VERSION).tar.gz >$@

src: phantom.tar.gz
	tar xf $<
	mv phantom-$(VERSION) src


# Building the code into a static library
FCFLAGS += $(STOPCOND_CFLAGS) $(FORSOCKETS_CFLAGS) $(OPENMP_FCFLAGS)
LDLIBS += -lm $(STOPCOND_LIBS) $(FORSOCKETS_LIBS)
LDFLAGS += $(OPENMP_FCFLAGS)

CODELIB = src/bin/libphantom-amuse.a

.PHONY: $(CODELIB)
$(CODELIB): | src
	$(MAKE) -C src -j 1 libphantom-amuse SETUP=amuse SYSTEM=gfortran OPENMP=yes


# Building the workers

INTERFACE_FCFLAGS = $(FCFLAGS) -Isrc/build -Lsrc/bin -lphantom-amuse 

phantom_worker.f90: interface.py
	amusifier --type=f90 interface.py PhantomInterface -o $@

phantom_worker.o: phantom_worker.f90
	$(MPIFC) -c -o $@ $(FCFLAGS) $<

phantom_worker: $(CODELIB) interface.o phantom_worker.o
	$(MPIFC) -o $@ $(LDFLAGS) $^ $(LDLIBS) $(INTERFACE_FCFLAGS)

interface.o: interface.f90
	$(MPIFC) -o $@ -c $(INTERFACE_FCFLAGS) $<

%.o: %.f90 $(CODELIB)
	$(MPIFC) -o $@ -c $(FCFLAGS) -Isrc/build $<


# Which packages contain which workers?
amuse-phantom_contains: phantom_worker


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
	-$(MAKE) -C src clean
	rm -rf *.o *.mod *worker*

.PHONY: distclean
distclean: clean
	-$(MAKE) -C src distclean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache

