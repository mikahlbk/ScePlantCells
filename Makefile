CC=g++ -std=c++11

CFLAGS=-c -Wall

all: program

program: main.o coord.o node.o cell.o tissue.o
		mkdir -p ./Animation ./DataOutput
		$(CC) main.o coord.o node.o cell.o tissue.o -o program

main.o: main.cpp
		$(CC) $(CFLAGS) main.cpp

coord.o: coord.cpp
		$(CC) $(CFLAGS) coord.cpp

node.o: node.cpp
		$(CC) $(CFLAGS) node.cpp

cell.o: cell.cpp
		$(CC) $(CFLAGS) cell.cpp

tissue.o: tissue.cpp
		$(CC) $(CFLAGS) tissue.cpp

clean:
		rm -rf *o program ./Animation ./DataOutput
