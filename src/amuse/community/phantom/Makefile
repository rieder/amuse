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

all: phantom_worker 

ifdef DOWNLOAD_CODES
$(CODEDIR)/Makefile:
	make -C . download
else
$(CODEDIR)/Makefile:
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "DOWNLOAD_CODES is not set. Phantom will not be downloaded and built."
	@echo "If you do want Phantom, set DOWNLOAD_CODES to 1."
	@echo "bash> export DOWNLOAD_CODES=1"
	@echo "csh> setenv DOWNLOAD_CODES 1"
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@make -s --no-print-directory -C . raise_error
endif

download:
	$(RM) -Rf src
	mkdir src
	$(DOWNLOAD_FROM_WEB)

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ worker_code worker_code.f90
	$(RM) ${CODELIB}
	make -C ${CODEDIR} clean

$(CODELIB): $(CODEDIR)/Makefile
	make -C ${CODEDIR} libphantom SETUP=amuse SYSTEM=gfortran OPENMP=yes
	ar crs ${CODELIB} ${CODEDIR}/build/*.o

worker_code.f90: interface.py
	$(CODE_GENERATOR) --type=f90 interface.py PhantomInterface -o $@

phantom_worker: worker_code.f90 $(CODELIB) $(OBJS)
	$(MPIFC) $(FCFLAGS) $(FS_FLAGS) $< $(OBJS) $(CODELIB) $(FS_LIBS) -fopenmp -o $@

%.o: %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<
