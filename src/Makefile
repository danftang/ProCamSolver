# My Makefile

MAIN = main.o
OBJS = FundamentalMatrix.o
#OBJS = staticMeasurementMatrix.o staticMotionMatrix.o staticStrumTriggsSolver.o
#CPPFLAGS = -O3 -frepo
LIBS =

all: $(OBJS) $(MAIN)
	g++ $(OPTFLAGS) $(CPPFLAGS) $(OBJS) $(MAIN) $(LIBS)

clean:
	rm $(OBJS) $(MAIN)
