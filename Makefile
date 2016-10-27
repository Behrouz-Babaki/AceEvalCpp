.PHONY: all clean test

CXX = g++
CFLAGS = --std=c++0x -g -I./AceEvalCpp
LFLAGS = -g

all: test/test

test/main.o: test/main.cpp AceEvalCpp/AceEvalCpp.hpp
	$(CXX) -c $(CFLAGS) $< -o ./test/main.o

test/test: test/main.o
	$(CXX) test/main.o -o test/test $(LFLAGS)

test: test/test
	@./test/test ./test/toy_bn.net.ac ./test/toy_bn.net.lmap

clean: 
	@rm -rf *~ ./test/*.o ./test/test ./test/*~ ./core
