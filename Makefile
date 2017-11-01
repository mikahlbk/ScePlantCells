CC=g++ -std=c++11

CFLAGS=-c -Wall

all: program

program: folder main.o coord.o node.o side.o cell.o cell_div.o tissue.o
		$(CC) main.o coord.o node.o side.o cell.o cell_div.o tissue.o -o program

folder: 
		mkdir -p ./Animation ./DataOutput

main.o: main.cpp
		$(CC) $(CFLAGS) main.cpp

coord.o: coord.cpp
		$(CC) $(CFLAGS) coord.cpp

node.o: node.cpp
		$(CC) $(CFLAGS) node.cpp

side.o: side.cpp
		$(CC) $(CFLAGS) side.cpp

cell.o: cell.cpp
		$(CC) $(CFLAGS) cell.cpp

cell_div.o: cell_div.cpp
		$(CC) $(CFLAGS) cell_div.cpp

tissue.o: tissue.cpp
		$(CC) $(CFLAGS) tissue.cpp

clean: wipe
		rm -rf cell_vec.txt side_vecs.txt  *o program

wipe:
		rm -rf ./Animation ./DataOutput
