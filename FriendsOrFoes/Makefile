CXX = g++
CXXLD = g++
CXXFLAGS = -O2
CXXCOMPILE = $(CXX) $(CXXFLAGS)
CXXLINK = $(CXXLD) $(CXXFLAGS) -o $@
INCLUDES = -I/usr/include/g++ -I../lib/poco/Util/include -I../include -I../include
LIB_PATH = -L../lib/poco/lib/Linux/x86_64 -L/usr/lib/sprng
LIBS = -lPocoFoundation -lPocoUtil -lm 

SOURCES = friends_or_foes.cpp population.cpp cell_cycle.cpp probability.cpp ../lib/RngStream/RngStream.cpp

OBJECTS = friends_or_foes.o population.o cell_cycle.o probability.o ../lib/RngStream/RngStream.o

all: friends_or_foes

friends_or_foes.o: friends_or_foes.cpp friends_or_foes.h clone.h

population.o: population.cpp population.h cell.h clone.h replication_record.h

cell_cycle.o: cell_cycle.cpp cell_cycle.h probability.h

probability.o: probability.h

friends_or_foes: $(OBJECTS)
	@rm -f friends_or_foes
	$(CXXLINK) $(OBJECTS) $(INCLUDE_PATH) $(DEFINES) $(LIB_PATH) $(LIBS)

.cpp.o:
	$(CXXCOMPILE) $(INCLUDES) $(DEFINES) -c -o $@ $<

clean:
	@rm -f friends_or_foes $(OBJECTS)

