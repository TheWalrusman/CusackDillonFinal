.PHONY: all clean

#macros
CXX = /usr/bin/g++-7
CXXFLAGS = -Wpedantic -Wall -Wextra -Wfloat-conversion -Werror
SOURCES = $(wildcard *.cpp)
HEADERS = $(wildcard *.h)
OBJECTS = $(SOURCES:%.cpp=%.o)

default: driver

#how to make a .o from a .cpp
%.o: %.cpp
	@$(CXX) $(CXXFLAGS) -c $< -o $@

#create an excutable called driver
driver: $(OBJECTS)
	@$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@

#remove core dumps, executables, depend file, and .o files
clean:
	-@rm -f core
	-@rm -f driver
	-@rm -f depend
	-@rm -f $(OBJECTS)

# Automatically generate dependencies and include them in Makefile
depend: $(SOURCES) $(HEADERS)
	@$(CXX) -MM *.cpp > $@

-include depend