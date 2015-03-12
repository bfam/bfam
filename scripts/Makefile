#
# This is an example Makefile to build a executable depending on bfam which
# uses bfam's exported Makefile.
#

include $(BFAM_DIR)/lib/bfam/bfam.mk

CC=$(BFAM_CC)
CFLAGS+=$(BFAM_INCLUDES) $(BFAM_CFLAGS)
LDFLAGS+=$(BFAM_LDFLAGS) $(BFAM_CFLAGS)
LDLIBS+=$(BFAM_LIBRARIES)

NDEBUG ?= 1
ifeq (NDEBUG, 1)
  CFLAGS+=$(BFAM_CFLAGS_DEBUG)
else
  CFLAGS+=$(BFAM_CFLAGS_RELEASE)
endif

ifeq ($(BFAM_HAVE_CUDA),ON)
  CFLAGS+=-I$(BFAM_CUDA_INCLUDE_DIRS)
endif
ifeq ($(BFAM_HAVE_OPENCL),ON)
  CFLAGS+=-I$(BFAM_OPENCL_INCLUDE_DIRS)
endif

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
  LDLIBS+=-lrt
endif

default: print_info main

print_info:
	@echo "\nFound BFAM with details: "
	@echo "   BFAM_VERSION             = $(BFAM_VERSION)"
	@echo "   BFAM_CC                  = $(BFAM_CC)"
	@echo "   BFAM_CXX                 = $(BFAM_CXX)"
	@echo "   BFAM_CFLAGS              = $(BFAM_CFLAGS)"
	@echo "   BFAM_CFLAGS_DEBUG        = $(BFAM_CFLAGS_DEBUG)"
	@echo "   BFAM_CFLAGS_RELEASE      = $(BFAM_CFLAGS_RELEASE)"
	@echo "   BFAM_INCLUDES            = $(BFAM_INCLUDES)"
	@echo "   BFAM_LIBRARIES           = $(BFAM_LIBRARIES)"
	@echo "   BFAM_LDFLAGS             = $(BFAM_LDFLAGS)"
	@echo "   BFAM_CUDA_INCLUDE_DIRS   = $(BFAM_CUDA_INCLUDE_DIRS)"
	@echo "   BFAM_OPENCL_INCLUDE_DIRS = $(BFAM_OPENCL_INCLUDE_DIRS)"
	@echo "   BFAM_HAVE_CUDA           = $(BFAM_HAVE_CUDA)"
	@echo "   BFAM_HAVE_OCCA           = $(BFAM_HAVE_OCCA)"
	@echo "   BFAM_HAVE_OPENCL         = $(BFAM_HAVE_OPENCL)"
	@echo "End of BFAM details\n"

main: main.o

.PHONY: clean
clean:
	rm -f main *.o
