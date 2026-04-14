AMUSE_VERSION := $(patsubst v%,%,$(shell git describe --tags))

ifeq (,$(AMUSE_VERSION))
    H := #
    AMUSE_VERSION := $(shell grep -v '^$H' ../../VERSION)
endif

export AMUSE_VERSION

