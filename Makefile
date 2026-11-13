PySRCS := $(wildcard *.py)

FLAKE8 = flake8 --max-line-length=100 --ignore=E221,E241,W503

OPT   = -g -O3 -Wall
CFLAGS = $(OPT)  -D_GNU_SOURCE -fPIC -fopenmp
LDFLAGS = -shared -fopenmp
LIBS = -lm
CC    = gcc

default:lib_radcouple.so check

lib_radcouple.so:lib_radcouple.c lib_vec.c lib_radcouple.h lib_vec.h
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ lib_radcouple.c lib_vec.c $(LIBS)

check:
	@for i in $(PySRCS); do echo  $(FLAKE8) $$i; $(FLAKE8) $$i; done;
