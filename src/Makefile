# My Makefile

MAIN = main.o
OBJS = transformMatrix.o FSolver.o PSolver.o SceneSynth.o RPSolver.o
CPPFLAGS = -O3 -I../
LIBS =

all: $(OBJS) $(MAIN)
	g++ $(OPTFLAG) $(CPPFLAGS) $(OBJS) $(MAIN) $(LIBS)

clean:
	rm $(OBJS) $(MAIN)
