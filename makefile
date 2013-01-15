CC = gcc
CFLAGS = -c -O2 -Wall
LDFLAGS = -lm -Wall -O2
FCC = gfortran
SOURCES1 = fastpopsim.c noia.c reg.c
SOURCES2 = constant_NSWC.f90 inc_beta.f90 fprob.f90
OBJECTS1 = $(SOURCES1:.c=.o)
UNAME := $(shell uname)
OBJECTS2 = $(SOURCES2:.f90=.o)
EXE1 = epiFit
SO = fprob.so

all: $(SOURCES1) $(SOURCES2) $(EXE1)

$(EXE1): $(OBJECTS1)
	$(CC) $(LDFLAGS) $(OBJECTS1) -ldl -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

ifeq ($(UNAME),Darwin)
all:
	$(FCC) -dynamiclib $(SOURCES2) -o $(SO)
	rm $(OBJECTS1)
else
all:
	$(FCC) -fPIC -c $(SOURCES2)
	$(FCC) -shared $(OBJECTS2) -o $(SO)
	rm $(OBJECTS1) $(OBJECTS2)
endif
clean:
	rm $(OBJECTS1) $(EXE1)

