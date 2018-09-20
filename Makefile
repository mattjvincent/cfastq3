TOP := $(shell /bin/pwd)

LIBS   = -lm -lz 
#CC     = c++
CC     = gcc
CFLAGS = -g -Wall
RM     = rm -f
INC    = -I$(TOP)


default: all

all: cfastq3

<<<<<<< HEAD
cfastq3: cfastq3.c
	$(CC) $(CFLAGS) $(LIBS) -o cfastq3 cfastq3.c
=======
cfastq3 : cfastq3.c
	$(CC) $(INC) $(CFLAGS) $(LIBS) libbloom/build/libbloom.a -o cfastq3 cfastq3.c
>>>>>>> master

clean :
	$(RM) cfastq3
	$(RM) -r cfastq3.dSYM

bloom : 
	(cd libbloom && make)
	

