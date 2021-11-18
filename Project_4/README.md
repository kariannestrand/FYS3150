## The Ising model

By using the Metropolis algorithm, we have implemented the Ising model in two dimensions two investigate phase transitions. The code is organized in two source files, functions.cpp and main.cpp and one header file, functions.hpp. 

In functions.hpp we have decleared all functions used in our code, and imported the libraries needed. functions.cpp contains all the calculations needed, and the Metropolis algrothim, while in main.cpp we have set all initial values. 

We have parallelized our code using OpenMP, at the temperature level, meaning that we parallelized the for loop over temperature. 

The lattice size L is an input argument to be given in the terminal. If the spins should start from ordered or unordered initial states can be chosen in line x in main.cpp. There is also the possibility to choose to run for one value of the temperature, or run for a vector of temperatures. If line x in main.cpp is set to true, the temperature is given as one value in the command line. If it is set to false, the code is running for a vector of temperatures that is set in main.cpp. The results are written to file if line x or x in main.cpp is set to true, depending on how many temperature values you want to run the code for. 

To build the code with a Linux terminal, write the following in the command line:

      g++ -O3 -larmadillo -std=c++11 *.cpp -fopenmp -o main.exe
  
To build the code with a MacOS terminal, write the following in the command line:

      g++ -O3 -std=c++11 *.cpp -Xpreprocessor -fopenmp -o main.exe -larmadillo -lomp
      
To run the code, in both cases, write

      ./main.exe
      
 in the command line, followed by your choise for the lattice size L and the temperature T if you run for one temperature. 
