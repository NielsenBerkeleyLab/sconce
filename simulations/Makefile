# compiler
CC=gcc
# compiler flags
CFLAGS=-lm


all: sconce_sim

sconce_sim: CNV.c
	$(CC) CNV.c ranlib.c rnglib.c -o sconce_sim $(CFLAGS)

.PHONY: clean
clean:
	rm -f sconce_sim

