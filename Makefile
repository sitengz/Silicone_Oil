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
MOLECULE ?= 1
PHASE ?= prod
BLOCK_NS ?= 5

.PHONY: all generate atsc4i surface-tension test clean

all: $(GENERATOR)

$(GENERATOR): $(SOURCES)
	$(CXX) $(CPPFLAGS) -std=c++17 $(CXXFLAGS) Generator/oil_generator.cpp $(LDFLAGS) $(LDLIBS) -o $@

generate: $(GENERATOR)
	./$(GENERATOR) --config "$(CONFIG)" $(GENERATOR_ARGS)

atsc4i:
	$(PYTHON) Analysis/atsc4i.py --config "$(CONFIG)" --molecule-id "$(MOLECULE)"

surface-tension:
	$(PYTHON) Analysis/surface_tension.py --config "$(CONFIG)" --phase "$(PHASE)" --block-ns "$(BLOCK_NS)"

test:
	$(PYTHON) -m unittest discover -s tests -v

clean:
	rm -f $(GENERATOR)
