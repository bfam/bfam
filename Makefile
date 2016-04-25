CC = mpicc
CXX = mpicxx
F77 = mpif77
FC = mpif90
CFLAGS = --std=gnu99 -g -O2

CFLAGS +=-fno-common -fomit-frame-pointer

CFLAGS += -Wconversion -Wno-sign-conversion \
          -Wcast-align -Wchar-subscripts -Wall -W \
          -Wpointer-arith -Wwrite-strings -Wformat-security -pedantic \
          -Wextra -Wno-unused-parameter -fPIC

LDFLAGS += -shared

# list of libraries to build
ALL_TPLS ?= third_party/p4est third_party/lua third_party/occa
TPLS ?= third_party/p4est

UNAME_S := $(shell uname -s)

# bfam flags

BFAM_HEADERS = src/bfam.h
BFAM_SOURCE = src/bfam.c
CPPFLAGS += -Isrc

# p4est flags
CPPFLAGS += -Ithird_party/p4est/local/include
LDFLAGS += -Lthird_party/p4est/local/lib
LDLIBS += -lp4est -lsc
ifeq ($(UNAME_S),Linux)
 LDFLAGS += -Wl,-rpath=$(CURDIR)/third_party/p4est/local/lib,--enable-new-dtags
endif

ifndef RELEASE
	CPPFLAGS += -DBFAM_DEBUG
	P4EST_DEBUG = --enable-debug
endif

ifdef USE_LUA
	CPPFLAGS += -DBFAM_USE_LUA

	TPLS += third_party/lua

	CPPFLAGS += -Ithird_party/lua/src
	LDFLAGS += -Lthird_party/lua/src
	LDLIBS += -llua -lm
endif

ifdef USE_BFAMO
	CPPFLAGS += -DBFAM_USE_BFAMO
	TPLS += third_party/occa
	CPPFLAGS += -Ithird_party/occa/include

	LDFLAGS += -Lthird_party/occa/lib
	LDLIBS += -locca

ifeq ($(UNAME_S),Linux)
	LDFLAGS += -Wl,-rpath=$(CURDIR)/third_party/occa/lib,--enable-new-dtags
endif

ifneq ($(BFAMO_REAL_DOUBLE), 0)
	CPPFLAGS += -DBFAMO_REAL_DOUBLE
endif

ifdef KERNEL_TIME
ifneq ($(KERNEL_TIME), 0)
	CPPFLAGS += -DKERNEL_TIME
endif
endif

ifdef COMPUTE_STATS
ifneq ($(COMPUTE_STATS), 0)
	CPPFLAGS += -DCOMPUTE_STATS
endif
endif

endif

SHAREDLIBRARY = libbfam2d.so libbfam3d.so

all:

tpls: $(ALL_TPLS)

third_party/lua:
	tar zxf third_party/lua-5.3.2.tar.gz
	mv lua-5.3.2 third_party/lua
	cd third_party/lua && $(MAKE) CC=$(CC) posix

third_party/occa:
	tar zxf third_party/occa-43903d8.tar.gz
	mv occa-master third_party/occa
	cd third_party/occa && $(MAKE) CC=$(CC) CXX=$(CXX) FC=$(FC)

third_party/p4est:
	tar zxf third_party/p4est-1.1.tar.gz
	mv p4est-1.1 third_party/p4est
	cd third_party/p4est \
		&& ./configure CC=$(CC) CXX=$(CXX) F77=$(F77) FC=$(FC) --enable-mpi \
		--without-blas $(P4EST_DEBUG) && $(MAKE) install

# Dependencies

shared: $(SHAREDLIBRARY)

bfam2d.o: CPPFLAGS+=-DBFAM_DGX_DIMENSION=2
bfam2d.o: $(BFAM_SOURCE) | $(TPLS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c $^ -o $@

libbfam2d.so: bfam2d.o
	$(CC) $(LDFLAGS) $(TARGET_ARCH) $(LOADLIBES) $(LDLIBS) $^ -o $@

bfam3d.o: CPPFLAGS+=-DBFAM_DGX_DIMENSION=3
bfam3d.o: $(BFAM_SOURCE) | $(TPLS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c $^ -o $@

libbfam3d.so: bfam3d.o
	$(CC) $(LDFLAGS) $(TARGET_ARCH) $(LOADLIBES) $(LDLIBS) $^ -o $@

# Rules
.PHONY: clean realclean
clean:
	rm -rf $(SHAREDLIBRARY) *.o

realclean: clean
	rm -rf $(ALL_TPLS)
# 	git clean -x -d -f

