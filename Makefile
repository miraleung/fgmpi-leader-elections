MPICC := mpicc -g -D_REENTRANT -W -Wall -O3
MPILD := $(MPICC) -lm

INC = 

CFILES := $(wildcard *.c)
APPSOBJS := $(patsubst %.c, %.o, $(CFILES)) 
HEADERS := $(wildcard *.h)

APPS = $(patsubst %.c, %, $(CFILES))

all: $(APPS)
	@echo Making $(APPS) ....

%.o: %.c $(HEADERS) Makefile
	$(MPICC) $(INC) -o $@ -c $<


$(APPS) : % : %.o $(MAKEDEPS) $(APPSOBJS) $(HEADERS)
	$(MPILD) $(INC) $(patsubst %, %.o, $@) -o $@ 



clean:
	rm -f *.o *.a core $(APPS)


FORCE:

first_target: all

.PHONY: all clean 
