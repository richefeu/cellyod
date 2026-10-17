# cellyod : C++17 sans dépendance ; l'interface (cellyod-gui) dépend en plus de raylib (brew install raylib)
CXX      ?= c++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra

CORE_SRC = src/params.cpp src/tissue.cpp src/io.cpp
CORE_HDR = src/geom.hpp src/params.hpp src/tissue.hpp src/io.hpp

RAYLIB = $(shell pkg-config --cflags --libs raylib)
ifeq ($(shell uname -s),Darwin)
	RAYLIB += -framework Cocoa -framework IOKit -framework OpenGL
endif

all: cellyod

cellyod: src/cli.cpp $(CORE_SRC) $(CORE_HDR)
	$(CXX) $(CXXFLAGS) -o $@ src/cli.cpp $(CORE_SRC)

clean:
	rm -f cellyod cellyod-gui

.PHONY: all clean
