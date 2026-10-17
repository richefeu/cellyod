# cellyod : C++17 sans dépendance ; l'interface (cellyod-gui) dépend en plus de raylib (brew install raylib)
CXX      ?= c++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra

CORE_SRC = src/params.cpp src/expr.cpp src/tissue.cpp src/io.cpp
CORE_HDR = src/geom.hpp src/expr.hpp src/params.hpp src/tissue.hpp src/io.hpp

RAYLIB = $(shell pkg-config --cflags --libs raylib)
ifeq ($(shell uname -s),Darwin)
	RAYLIB += -framework Cocoa -framework IOKit -framework OpenGL
endif

all: cellyod cellyod-gui

cellyod: src/cli.cpp $(CORE_SRC) $(CORE_HDR)
	$(CXX) $(CXXFLAGS) -o $@ src/cli.cpp $(CORE_SRC)

# raygui.h et rini.h (single-header) déclenchent des avertissements : neutralisés pour cette cible seulement
cellyod-gui: gui/cellyod_gui.cpp $(CORE_SRC) $(CORE_HDR) gui/raygui.h gui/rini.h gui/guifont_ttf.h
	$(CXX) $(CXXFLAGS) -Wno-missing-field-initializers -Wno-unused-parameter -Wno-unused-function \
		-Wno-enum-compare -Wno-deprecated-declarations -Wno-unused-but-set-variable -o $@ gui/cellyod_gui.cpp $(CORE_SRC) $(RAYLIB)

clean:
	rm -f cellyod cellyod-gui

.PHONY: all clean
