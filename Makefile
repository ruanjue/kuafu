VERSION=1.0

CC := gcc

ifeq (0, ${MAKELEVEL})
TIMESTAMP=$(shell date)
endif

ifeq (1, ${DEBUG})
CFLAGS=-g3 -ggdb -W -Wall -Wno-unused-but-set-variable -O0 -DTIMESTAMP="$(TIMESTAMP)" -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE -mpopcnt -msse4.2
else
CFLAGS=-g3 -W -Wall -Wno-unused-but-set-variable -O4 -DTIMESTAMP="$(TIMESTAMP)" -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE -mpopcnt -msse4.2
endif

INSTALLDIR=/usr/local/bin
GLIBS=-lm -lrt -lpthread -lz
GENERIC_SRC=mem_share.h chararray.h pgzf.h filereader.h bitvec.h hashset.h sort.h list.h dna.h thread.h

PROGS=kuafu-sra kuafu-pog

all: $(PROGS)

kuafu-sra: $(GENERIC_SRC) srb.h sra.h kuafu-sra.c
	$(CC) $(CFLAGS) -o $@ kuafu-sra.c $(GLIBS)

kuafu-pog: $(GENERIC_SRC) srb.h pog.h kuafu-pog.c
	$(CC) $(CFLAGS) -o $@ kuafu-pog.c $(GLIBS)

clean:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out $(PROGS)

clear:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out

install:
	cp -fvu $(PROGS) $(INSTALLDIR)
