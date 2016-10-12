TARGET=./test_driver
SRC=$(wildcard *.c) perf/perf.c
CFLAGS+=-Wall -Wextra -std=gnu99
LDFLAGS+=

ifdef DEBUG
CFLAGS+=-O0 -ggdb #-fsanitize=address -fsanitize=undefined
LDFLAGS+= #-fsanitize=address -fsanitize=undefined
else
CFLAGS+=-O3
LDFLAGS+=
endif


ifdef USE_MKL
CFLAGS+=-DUSE_MKL=1
EXTRA_OBJ=${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a
LDFLAGS+= -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
	${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a \
	-Wl,--end-group ${MKLROOT}/lib/intel64/libmkl_blacs_openmpi_lp64.a \
	-ldl -lpthread -lm -fopenmp
endif


OBJ=$(SRC:.c=.o)
DEP=$(SRC:.c=.d)

all: $(TARGET)

-include $(DEP)

$(TARGET): $(OBJ)
	$(CC) $^ -o $@ $(LDFLAGS)

%.o: %.c
	@$(CC) -MM $(CFLAGS) $*.c > $*.d
	$(CC) -c $(CFLAGS) $*.c -o $*.o

clean:
	$(RM) $(OBJ) $(DEP) *.d *.o

mrproper: clean
	$(RM) $(TARGET)

test:
	$(TARGET)
