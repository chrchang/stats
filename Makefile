CWARN = -Wall -Wextra -Wshadow -Wformat-security -Wstrict-aliasing -Wdouble-promotion -Wfloat-conversion
CXXWARN = ${CWARN} -Wold-style-cast

ISRC = include/plink2_base.cc \
       include/plink2_float.cc \
       include/plink2_highprec.cc \
       include/plink2_ln.cc

IHDR = $(ISRC:.cc=.h)

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

%.o: %.cc $(IHDR)
	g++ -c $(CXXFLAGS) -o $@ $<

all: fisher_test hwe_test

fisher_test: $(OBJ) include/fisher.o fisher.o fisher_test.o
	g++ $(OBJ) include/fisher.o fisher.o fisher_test.o -o fisher_test $(LINKFLAGS)

hwe_test: $(OBJ) include/plink2_hwe.o hwe_test.o
	g++ $(OBJ) include/plink2_hwe.o hwe_test.o -o hwe_test $(LINKFLAGS)

.PHONY: clean
clean:
	rm -f *.o
	rm -f include/*.o
	rm -f fisher_test
	rm -f hwe_test
