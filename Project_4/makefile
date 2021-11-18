all: build run

build:
	g++ -O3 -larmadillo -std=c++11 *.cpp -fopenmp -o main.exe

run:
	./main.exe

macos:
	g++ -O3 -std=c++11 *.cpp -Xpreprocessor -fopenmp -o main.exe -larmadillo -lomp
