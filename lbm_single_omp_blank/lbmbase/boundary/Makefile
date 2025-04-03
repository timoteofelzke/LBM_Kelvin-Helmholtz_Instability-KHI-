CPPS = $(wildcard *.cpp)
OBJS = $(CPPS:.cpp=.o)

all : $(OBJS)

%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	\rm *.o *~ p1 *.d -f
