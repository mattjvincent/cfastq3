TOP := $(shell /bin/pwd)

LIBS   = -lm -lz 
#CC     = c++
CC     = gcc
CFLAGS = -g -Wall
RM     = rm -f
INC    = -I$(TOP)


default: all

all: cfastq3

cfastq3: cfastq3.c
	$(CC) $(CFLAGS) $(LIBS) -o cfastq3 cfastq3.c

clean :
	$(RM) cfastq3
	$(RM) -r cfastq3.dSYM


