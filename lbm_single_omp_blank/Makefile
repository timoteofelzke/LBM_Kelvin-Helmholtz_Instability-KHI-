export CC=g++
export CFLAGS=-Wall -O2 -msse4.2 -fomit-frame-pointer -ffast-math -funroll-all-loops -mfpmath=sse -finline-limit=800 -Wno-unknown-pragmas -MMD -fopenmp
export LFLAGS=-lboost_program_options -lboost_timer -lboost_system -lz -fopenmp

SUBDIRS      := lbmbase problems
SRC          := $(foreach sdir,$(SUBDIRS),$(wildcard $(sdir)/*.cpp))
NOBJ	     := $(foreach sdir,$(SUBDIRS),$(sdir)/$(sdir).o)
OBJ          := $(SRC:.cpp=.o)
LOCALCPPS    := $(wildcard *.cpp)
LOCALOBJS    := $(LOCALCPPS:.cpp=.o)
LOCALDEPS    := $(LOCALOBJS:.o=.d)

EXECUTABLE=grad-lbm

all: $(LOCALOBJS) $(SUBDIRS)
	$(CC) $(NOBJ) $(LOCALOBJS) $(LFLAGS) -o $(EXECUTABLE)

-include $(LOCALDEPS)

%.o : %.cpp %.d
	$(CC) $(CFLAGS) -c $< -o $@

%.d : %.cpp
	$(CC) -E -MMD $< >> /dev/null

$(SUBDIRS):
	$(MAKE) -C $@

debug : CFLAGS = -DDEBUG -g -O0 -gdwarf-2 -MMD
debug : LFLAGS += -DDEBUG -g -gdwarf-2
debug : all

clean:
	\rm *.o *.d *~ grad-lbm p1 -f
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean; \
	done 

.PHONY: all $(SUBDIRS)
