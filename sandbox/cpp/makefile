# Variables 
CC = g++
CFLAGS = -Wall -g -std=c++11

# Executable
main: main.o genSpectra.o 
	$(CC) $(CFLAGS) -o main main.o genSpectra.o 

# Object files
main.o: main.cpp genSpectra.hpp
	$(CC) $(CFLAGS) -c main.cpp loadFasta.cpp

genSpectra.o: genSpectra.hpp 
	$(CC) $(CFLAGS) -c genSpectra.cpp 

clean:
	rm main *.o