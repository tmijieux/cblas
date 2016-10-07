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
