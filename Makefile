CXX ?= g++
PYTHON ?= python3
CXXFLAGS ?= -O2 -Wall -Wextra -Wpedantic
CPPFLAGS ?=
LDFLAGS ?=
LDLIBS ?=

GENERATOR := Generator/oil_generator
SOURCES := Generator/oil_generator.cpp Generator/config_input.hpp
CONFIG ?= model.conf
GENERATOR_ARGS ?=

.PHONY: all generate test clean

all: $(GENERATOR)

$(GENERATOR): $(SOURCES)
	$(CXX) $(CPPFLAGS) -std=c++17 $(CXXFLAGS) Generator/oil_generator.cpp $(LDFLAGS) $(LDLIBS) -o $@

generate: $(GENERATOR)
	./$(GENERATOR) --config "$(CONFIG)" $(GENERATOR_ARGS)

test:
	$(PYTHON) -m unittest discover -s tests -v

clean:
	rm -f $(GENERATOR)
