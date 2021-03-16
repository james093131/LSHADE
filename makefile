all:
	g++ -Wall -O3 -o LSHADE main.cpp

clean:
	rm -f main *.o
