# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../../../../../../Env/Amuse/share/amuse
-include $(AMUSE_DIR)/config.mk

MPICXX   ?= mpicxx

CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODELIB = src/libgasoline.a

CODE_GENERATOR = $(AMUSE_DIR)/build.py

DOWNLOAD_FROM_WEB = $(PYTHON) ./download_http.py


all: gasoline_worker 

ifdef DOWNLOAD_CODES
$(CODEDIR)/Makefile:
	make -c . download
else
$(CODEDIR)/Makefile:
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@echo "DOWNLOAD_CODES is not set. Gasoline will not be downloaded and built."
	@echo "If you do want Gasoline, set DOWNLOAD_CODES to 1."
	@echo "bash> export DOWNLOAD_CODES=1"
	@echo "csh> setenv DOWNLOAD_CODES 1"
	@echo ""
	@echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	@echo ""
	@make -s --no-print-directory -C . raise_error
endif

download:
	$(RM) -Rf .pc
	$(RM) -Rf src
	mkdir src
	$(DOWNLOAD_FROM_WEB)
	#$(PATCH_FILES)

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ gasoline_worker worker_code.cc
	$(RM) -Rf $(CODEDIR)/src/*.o
	$(RM) -Rf $(CODEDIR)/src/*.so
	$(RM) -Rf src/libgasoline.a
	make -C src clean

$(CODELIB):
	make -C src all

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py GasolineInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py GasolineInterface -o $@

gasoline_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $< $(OBJS) $(CODELIB) -o $@

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< 
