# My Makefile

MAIN = main.o
OBJS = src/FundamentalMatrix.o
#CPPFLAGS = -O3 -frepo
LIBS =

all: $(OBJS) $(MAIN)
	g++ $(OPTFLAGS) $(CPPFLAGS) $(OBJS) $(MAIN) $(LIBS)

clean:
	rm $(OBJS) $(MAIN)
