mac: compile_m link_m execute
linux: compile_l link_l execute

compile_m:
	g++ -c *.cpp -std=c++11

link_m:
	g++ -o main.exe *.o -larmadillo


compile_l:
	g++ -c *.cpp

link_l:
	g++ -o main.exe *.o

execute:
	./main.exe
