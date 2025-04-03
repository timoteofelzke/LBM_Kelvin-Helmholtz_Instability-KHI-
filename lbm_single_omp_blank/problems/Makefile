DIRNAME=problems
CPPS = $(wildcard *.cpp)
OBJS = $(CPPS:.cpp=.o)
DEPS = $(OBJS:.o=.d)

all : $(OBJS)
	ld -r $(OBJS) -o $(DIRNAME).o

-include $(DEPS)

%.o : %.cpp %.d
	$(CC) $(CFLAGS) -c $< -o $@

%.d : %.cpp
	$(CC) -E -MMD $< >> /dev/null

clean:
	\rm *.o *~ *.d p1 -f
