CXX = g++
CXXLD = g++
CXXFLAGS = -O2
CXXCOMPILE = $(CXX) $(CXXFLAGS)
CXXLINK = $(CXXLD) $(CXXFLAGS) -o $@
INCLUDES = -I/usr/include/g++ -I../lib/poco/Util/include -I../include -I../include
LIB_PATH = -L../lib/poco/lib/Linux/x86_64
LIBS = -lPocoFoundation -lPocoUtil -lm 

COMMON_SOURCES = friends_or_foes.cpp population.cpp

RUN_SOURCES = probability.cpp cell_cycle.cpp ../lib/RngStream/RngStream.cpp run_friends_or_foes.cpp

COMMON_OBJECTS = friends_or_foes.o population.o
RUN_OBJECTS = probability.o cell_cycle.o ../lib/RngStream/RngStream.o run_friends_or_foes.o

all: run_friends_or_foes

run_friends_or_foes.o: run_friends_or_foes.cpp run_friends_or_foes.h

run_friends_or_foes.h: friends_or_foes.h

friends_or_foes.o: friends_or_foes.cpp friends_or_foes.h clone.h

population.o: population.cpp population.h cell.h clone.h replication_record.h

cell_cycle.o: cell_cycle.cpp cell_cycle.h probability.h

probability.o: probability.h

run_friends_or_foes: $(COMMON_OBJECTS) $(RUN_OBJECTS)
	@rm -f run_friends_or_foes
	$(CXXLINK) $(COMMON_OBJECTS) $(RUN_OBJECTS) $(INCLUDE_PATH) $(DEFINES) $(LIB_PATH) $(LIBS)

.cpp.o:
	$(CXXCOMPILE) $(INCLUDES) $(DEFINES) -c -o $@ $<

clean:
	@rm -f run_friends_or_foes $(OBJECTS)

