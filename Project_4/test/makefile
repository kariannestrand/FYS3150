all: compile link run

compile:
	g++ -c *.cpp -std=c++11 -O3

link:
	g++ -o main.exe *.o -larmadillo

run:
	./main.exe