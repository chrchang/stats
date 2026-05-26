CWARN = -Wall -Wextra -Wshadow -Wformat-security -Wstrict-aliasing -Wdouble-promotion -Wfloat-conversion
CXXWARN = ${CWARN} -Wold-style-cast

ISRC = include/binom.cc \
       include/fisher.cc \
       include/plink2_base.cc \
       include/plink2_float.cc \
       include/plink2_highprec.cc \
       include/plink2_hwe.cc \
       include/plink2_ln.cc

CCHDR = $(ISRC:.cc=.h)

OBJ = $(ISRC:.cc=.o)

CINCLUDE =
CXXINCLUDE =

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

%.o: %.cc $(CCHDR)
	g++ -c $(CXXFLAGS) -o $@ $<

all: binom_demo fisher_demo hwe_demo

binom_demo: $(OBJ) binom_demo.o
	g++ $(OBJ) binom_demo.o -o binom_demo $(LINKFLAGS)

fisher_demo: $(OBJ) fisher_demo.o
	g++ $(OBJ) fisher_demo.o -o fisher_demo $(LINKFLAGS)

hwe_demo: $(OBJ) hwe_demo.o
	g++ $(OBJ) hwe_demo.o -o hwe_demo $(LINKFLAGS)

.PHONY: clean
clean:
	rm -f *.o
	rm -f include/*.o
	rm -f binom_demo
	rm -f fisher_demo
	rm -f hwe_demo
