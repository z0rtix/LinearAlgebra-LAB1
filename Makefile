CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2
LDFLAGS =

SRCS = main.cpp matrix.cpp gauss.cpp LU.cpp tests.cpp
OBJS = $(SRCS:.cpp=.o)

TARGET = main

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean