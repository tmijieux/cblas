TARGET=test
SRC=$(wildcard *.c) perf/perf.c
CFLAGS=-Wall -Wextra -std=gnu99 
LDFLAGS=
ifdef DEBUG
CFLAGS+=-O0 -ggdb #-fsanitize=address -fsanitize=undefined
LDFLAGS+= #-fsanitize=address -fsanitize=undefined
else
CFLAGS+=-O3
endif

OBJ=$(SRC:.c=.o)

$(TARGET): $(OBJ)
	$(CC) $^ -o $@ $(LDFLAGS)


clean:
	$(RM) $(OBJ)
