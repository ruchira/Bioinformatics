CXX = mpicxx
CXXLD = mpicxx
CXXFLAGS = -g -DUSING_MPI
CXXCOMPILE = $(CXX) $(CXXFLAGS)
CXXLINK = $(CXXLD) $(CXXFLAGS) -o $@
INCLUDES = -I/usr/include/g++ -I../lib/poco/Util/include -I../include -I../lib/allegro/include -I../lib/allegro/build/include
LIB_PATH = -L../lib/poco/lib/Linux/x86_64 -L../lib/allegro/build/lib
LIBS = -lPocoFoundation -lPocoUtil -lm -lgsl -lgslcblas -lprotobuf
SEE_LIBS = -lalleg

COMMON_SOURCES = friends_or_foes.cpp hex_const.cpp hex_cell.cpp hex_replication_record.cpp population.cpp hex_population.cpp probability.cpp ../lib/RngStream/RngStream.cpp hex_friends_or_foes_history.proto

RUN_SOURCES = run_friends_or_foes.cpp cell_cycle.cpp hex_cell_cycle.cpp

SEE_SOURCES = see_friends_or_foes.cpp visualize_hex_population.cpp hexagon_renderings.cpp hexagon_rendering.cpp hsv_to_rgb.cpp distinct_hues.cpp

COMMON_OBJECTS = friends_or_foes.o hex_const.o hex_cell.o hex_replication_record.o population.o hex_population.o probability.o ../lib/RngStream/RngStream.o hex_friends_or_foes_history.o

RUN_OBJECTS = run_friends_or_foes.o cell_cycle.o hex_cell_cycle.o

SEE_OBJECTS = see_friends_or_foes.o visualize_hex_population.o hexagon_renderings.o hexagon_rendering.o hsv_to_rgb.o distinct_hues.o

all: run_friends_or_foes see_friends_or_foes

run_friends_or_foes.o: run_friends_or_foes.cpp run_friends_or_foes.h hex_replication_record.h

run_friends_or_foes.h: friends_or_foes.h cell_cycle.h

friends_or_foes.h: clone.h population.h

friends_or_foes.o: friends_or_foes.cpp friends_or_foes.h hex_population.h

hex_const.o: hex_const.cpp hex_const.h

hex_cell.o: hex_cell.cpp hex_cell.h hex_const.h

hex_cell.h: cell.h

cell.h: clone.h

hex_replication_record.o: hex_replication_record.cpp hex_replication_record.h hex_cell.h

hex_replication_record.h: replication_record.h

replication_record.h: cell.h

population.o: population.cpp population.h 

population.h: replication_record.h

hex_population.o: hex_population.cpp hex_population.h probability.h

hex_population.h: population.h hex_cell.h

probability.o: probability.cpp probability.h

cell_cycle.o: cell_cycle.cpp cell_cycle.h 

cell_cycle.h: population.h probability.h

hex_cell_cycle.o: hex_cell_cycle.cpp hex_cell_cycle.h

hex_cell_cycle.h: cell_cycle.h

hex_friends_or_foes_history.o: hex_friends_or_foes_history.proto
	protoc --cpp_out=. hex_friends_or_foes_history.proto; $(CXXCOMPILE) $(INCLUDES) $(DEFINES) -c hex_friends_or_foes_history.pb.cc -o hex_friends_or_foes_history.o

hex_friends_or_foes_history.pb.h: hex_friends_or_foes_history.proto
	protoc --cpp_out=. hex_friends_or_foes_history.proto

hsv_ro_rgb.o: hsv_to_rgb.cpp hsv_to_rgb.h 

see_friends_or_foes.h: friends_or_foes.h hex_friends_or_foes_history.pb.h hex_cell.h hex_replication_record.h

see_friends_or_foes.o: see_friends_or_foes.cpp see_friends_or_foes.h visualize_hex_population.h hexagon_renderings.h

visualize_hex_population.h: hex_population.h hexagon_rendering.h

visualize_hex_population.o: visualize_hex_population.cpp visualize_hex_population.h hsv_to_rgb.h distinct_hues.h

hexagon_renderings.h: hexagon_rendering.h

hexagon_renderings.o: hexagon_renderings.cpp hexagon_renderings.h

hexagon_rendering.o: hexagon_rendering.cpp hexagon_rendering.h

hsv_to_rgb.o: hsv_to_rgb.cpp hsv_to_rgb.h

distinct_hues.o: distinct_hues.cpp distinct_hues.h

run_friends_or_foes: $(COMMON_OBJECTS) $(RUN_OBJECTS)
	@rm -f run_friends_or_foes
	$(CXXLINK) $(COMMON_OBJECTS) $(RUN_OBJECTS) $(INCLUDE_PATH) $(DEFINES) $(LIB_PATH) $(LIBS)

see_friends_or_foes: $(COMMON_OBJECTS) $(SEE_OBJECTS)
	@rm -f see_friends_or_foes
	$(CXXLINK) $(COMMON_OBJECTS) $(SEE_OBJECTS) $(INCLUDE_PATH) $(DEFINES) $(LIB_PATH) $(LIBS) $(SEE_LIBS)

.cpp.o:
	$(CXXCOMPILE) $(INCLUDES) $(DEFINES) -c -o $@ $<

clean:
	@rm -f run_friends_or_foes see_friends_or_foes $(COMMON_OBJECTS) $(RUN_OBJECTS) $(SEE_OBJECTS) hex_cell_history.pb.h hex_cell_history.pb.cc

