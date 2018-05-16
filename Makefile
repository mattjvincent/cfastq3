LIBS   = -lm -lz
CC     = gcc
CFLAGS = -g -Wall
RM     = rm -f

default: all

all: cfastq3

cfastq3: cfastq3.c
	$(CC) $(CFLAGS) $(LIBS) -o test test.c
	$(CC) $(CFLAGS) $(LIBS) -o cfastq3 cfastq3.c

clean veryclean:
	$(RM) cfastq3
	$(RM) -r cfastq3.dSYM

