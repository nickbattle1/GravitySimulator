CXX      := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -Iinclude
LDFLAGS  := -lsfml-graphics -lsfml-window -lsfml-system

SOURCES := src/main.cpp src/CelestialBody.cpp src/GravitySimulator.cpp src/Scenarios.cpp
TARGET  := gravitysim

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f $(TARGET)

.PHONY: all clean