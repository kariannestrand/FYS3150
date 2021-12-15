## The Double slit experiment

We have used the Crank Nicolson method to solve the Schrödinger equation in two spatial and one time dimension. This was done to try to simulate results from the Double slit experiment. 

The code is organized in two source files, [functions.cpp](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_5/functions.cpp) and [main.cpp](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_5/main.cpp)  and one header file, [functions.hpp](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_5/functions.hpp) . 

In [functions.hpp](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_5/functions.hpp) we have decleared all functions used in our code, and imported the libraries needed. [functions.cpp](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_5/functions.cpp) contains all the calculations needed, while in [main.cpp](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_5/main.cpp) we have set all initial values. 

To build the code with a Linux terminal, write the following in the command line:

      g++ -O3 -larmadillo -std=c++11 *.cpp -o main.exe
  
To build the code with a MacOS terminal, write the following in the command line:

      g++ -O3 -std=c++11 *.cpp -o main.exe -larmadillo
      
To run the code, in both cases, write

      ./main.exe
      
 
