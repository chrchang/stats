CWARN = -Wall -Wextra -Wshadow -Wformat-security -Wstrict-aliasing -Wdouble-promotion -Wfloat-conversion
CXXWARN = ${CWARN} -Wold-style-cast

ISRC = include/plink2_base.cc \
       include/plink2_float.cc \
       include/plink2_highprec.cc \
       include/plink2_hwe.cc \
       include/plink2_ln.cc

GCSRC = mini-gmp/mini-gmp.c

OBJ = $(ISRC:.cc=.o) $(GCSRC:.c=.o)

CINCLUDE = -Imini-gmp
CXXINCLUDE = -Imini-gmp

CLEAN = *.o \
        include/*.o \
        mini-gmp/*.o

BASEFLAGS=-ffp-contract=off
CFLAGS=-O2 -std=gnu99 ${BASEFLAGS}
CXXFLAGS=-O2 -std=c++17 ${BASEFLAGS}
LINKFLAGS=

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
  CXXFLAGS=-O2 -std=c++17 -stdlib=libc++ ${BASEFLAGS}
endif

%.o: %.c
	gcc -c $(CFLAGS) -o $@ $<

%.o: %.cc
	g++ -c $(CXXFLAGS) -o $@ $<

all: hwe_test

hwe_test: $(OBJ) hwe_test.o
	g++ $(OBJ) hwe_test.o -o hwe_test $(LINKFLAGS)

.PHONY: clean
clean:
	rm -f *.o
	rm -f include/*.o
	rm -f hwe_test
