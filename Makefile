TARGET=driver pcache
driver_SRC= \
	perf/perf.c \
	daxpy.c \
	ddot.c \
	dgemm.c \
	dgemv.c \
	driver.c \
	util.c

pcache_SRC = \
	pcache.c \
	util.c 

CFLAGS+=-Wall -Wextra -std=gnu99 -fopenmp -march=native 
LDFLAGS+=-lm

ifdef DEBUG
CFLAGS+=-O0 -ggdb -fsanitize=address -fsanitize=undefined
LDFLAGS+= -fsanitize=address -fsanitize=undefined
else
CFLAGS+=-O3 -funroll-loops 
LDFLAGS+=
endif

ifdef USE_MKL
CFLAGS+=-DUSE_MKL=1
LDFLAGS+= -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
	${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a \
	-Wl,--end-group ${MKLROOT}/lib/intel64/libmkl_blacs_openmpi_lp64.a \
	-ldl -lpthread -lm -fopenmp
endif

all: $(TARGET)

-include $(DEP)

driver: $(driver_SRC) Makefile
	$(CC) $($@_SRC) -o $@ $(CFLAGS) $(LDFLAGS)

pcache: $(pcache_SRC) Makefile
	$(CC) $($@_SRC) -o $@ $(CFLAGS) $(LDFLAGS)


# %.o: %.c
# 	@$(CC) -MM $(CFLAGS) $*.c > $*.d
# 	$(CC) -c $(CFLAGS) $*.c -o $*.o

clean:
	$(RM) $(TARGET) $(OBJ) $(DEP) *.d *.o

mrproper: clean
	$(RM) $(TARGET)

test:
	$(TARGET)
