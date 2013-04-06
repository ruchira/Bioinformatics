CXX = mpicxx
CXXLD = mpicxx
CXXFLAGS = -O2 -DUSING_MPI
CXXCOMPILE = $(CXX) $(CXXFLAGS)
CXXLINK = $(CXXLD) $(CXXFLAGS) -o $@
INCLUDES = -I/usr/include/g++ -I../lib/poco/Util/include -I../include -I../include
LIB_PATH = -L../lib/poco/lib/Linux/x86_64
LIBS = -lPocoFoundation -lPocoUtil -lm 

COMMON_SOURCES = friends_or_foes.cpp population.cpp hex_const.cpp hex_cell.cpp hex_replication_record.cpp

RUN_SOURCES = probability.cpp cell_cycle.cpp ../lib/RngStream/RngStream.cpp run_friends_or_foes.cpp

COMMON_OBJECTS = friends_or_foes.o population.o hex_const.o hex_cell.o hex_replication_record.o
RUN_OBJECTS = probability.o cell_cycle.o ../lib/RngStream/RngStream.o run_friends_or_foes.o

all: run_friends_or_foes

run_friends_or_foes.o: run_friends_or_foes.cpp run_friends_or_foes.h

run_friends_or_foes.h: friends_or_foes.h

friends_or_foes.o: friends_or_foes.cpp friends_or_foes.h clone.h

hex_const.cpp: hex_const.h

hex_cell.cpp: hex_cell.h hex_const.h

hex_replication_record.h: replication_record.h

hex_replication_record.o: hex_replication_record.cpp hex_replication_record.h hex_cell.h

population.o: population.cpp population.h cell.h clone.h replication_record.h

cell_cycle.o: cell_cycle.cpp cell_cycle.h probability.h

probability.o: probability.h

run_friends_or_foes: $(COMMON_OBJECTS) $(RUN_OBJECTS)
	@rm -f run_friends_or_foes
	$(CXXLINK) $(COMMON_OBJECTS) $(RUN_OBJECTS) $(INCLUDE_PATH) $(DEFINES) $(LIB_PATH) $(LIBS)

.cpp.o:
	$(CXXCOMPILE) $(INCLUDES) $(DEFINES) -c -o $@ $<

clean:
	@rm -f run_friends_or_foes $(COMMON_OBJECTS) $(RUN_OBJECTS)

