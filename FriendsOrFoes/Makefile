CXX = g++
CXXLD = g++
CXXFLAGS = -O2
CXXCOMPILE = $(CXX) $(CXXFLAGS)
CXXLINK = $(CXXLD) $(CXXFLAGS) -o $@
INCLUDES = -I/usr/include/g++ -I../lib/poco/Util/include
LIB_PATH = -L../lib/poco/lib/Linux/x86_64
LIBS = -lPocoFoundation -lPocoUtil -lm

SOURCES = friends_or_foes.cpp

OBJECTS = friends_or_foes.o

all: friends_or_foes

friends_or_foes.cpp: friends_or_foes.h clone.h

friends_or_foes: $(OBJECTS)
	@rm -f friends_or_foes
	$(CXXLINK) $(OBJECTS) $(INCLUDE_PATH) $(DEFINES) $(LIB_PATH) $(LIBS)

.cpp.o:
	$(CXXCOMPILE) $(INCLUDE_PATH) $(DEFINES) -c -o $@ $<

clean:
	@rm -f friends_or_foes $(OBJECTS)
