# Compiler settings
CXX := g++
CXXFLAGS := -std=c++11 -Wall -Wextra

# Executable name
TARGET := afterglow

# Source file
SRCS := afterglow.cpp

# Object files
OBJS := $(SRCS:.cpp=.o)

# Build rule
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)

# Object file rule
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Plot rule
plot:
	python3 plotter2.py

# Clean rule
clean:
	rm -f $(OBJS) $(TARGET)
