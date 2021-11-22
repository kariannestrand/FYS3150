## The Ising model

By using the Metropolis algorithm, we have implemented the Ising model in two dimensions two investigate phase transitions. The code is organized in two source files, [functions.cpp](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_4/functions.cpp) and [main.cpp](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_4/main.cpp)  and one header file, [functions.hpp](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_4/functions.hpp) . 

In [functions.hpp](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_4/functions.hpp) we have decleared all functions used in our code, and imported the libraries needed. [functions.cpp](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_4/functions.cpp) contains all the calculations needed, and the Metropolis algrothim, while in [main.cpp](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_4/main.cpp) we have set all initial values. 

We have parallelized our code using OpenMP, at the temperature level, meaning that we parallelized the loop over temperature. We are also optimizing by using the -O3 flag.

The lattice size L is an input argument to be given in the terminal. If the spins should start from ordered or unordered initial states can be chosen in line 8 in main.cpp. There is also the possibility to choose to run for one value of the temperature, or run for a vector of temperatures. If line 9 in main.cpp is set to true, the temperature is given as one value in the command line. If it is set to false, the code is running for a vector of temperatures that is set in main.cpp. The results are written to file if line 10 in main.cpp is set to true.

To build the code with a Linux terminal, write the following in the command line:

      g++ -O3 -larmadillo -std=c++11 *.cpp -fopenmp -o main.exe
  
To build the code with a MacOS terminal, write the following in the command line:

      g++ -O3 -std=c++11 *.cpp -Xpreprocessor -fopenmp -o main.exe -larmadillo -lomp
      
To run the code, in both cases, write

      ./main.exe
      
 in the command line, followed by your choise for the lattice size L and the temperature T if you run for one temperature. 
