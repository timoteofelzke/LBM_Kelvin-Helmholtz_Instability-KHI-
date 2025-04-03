SUBDIRS      := lattice collision streaming vtkio force boundary
SRC          := $(foreach sdir,$(SUBDIRS),$(wildcard $(sdir)/*.cpp))
OBJ          := $(SRC:.cpp=.o)
DEP 	     := $(OBJ:.o=.d)
LOCALCPPS    := $(wildcard *.cpp)
LOCALOBJS    := $(LOCALCPPS:.cpp=.o)
LOCALDEPS    := $(LOCALOBJS:.o=.d)

all: $(LOCALOBJS) $(SUBDIRS)
	 ld -r $(OBJ) $(LOCALOBJS) -o lbmbase.o

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
	\rm *.o *.d *~ p1 -f
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean; \
	done 

.PHONY: all $(SUBDIRS)
