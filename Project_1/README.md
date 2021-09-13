# Numerical solution of the one-dimensional Poisson equation

We have organized the code for the project in two source files and one header file, [functions.cpp](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_1/functions.cpp), [functions.hpp](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_1/functions.hpp) and 
[main.cpp](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_1/main.cpp). 
[functions.hpp](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_1/functions.hpp) declares all the functions used for the project. 
[functions.cpp](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_1/functions.cpp) contains all the calculations. [main.cpp](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_1/main.cpp) defines all the initial values, does a timing test of the algorithms, writes the necessary output to file and calls all necessary functions.

To get the output written to file, line 46 in [main.cpp](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_1/main.cpp) must be set to "true" like displayed here:

    // Creates files exact_n.txt and approx_n.txt if true
    bool write_to_file = true;

If you do not want the output written to file, simply change "true" to "false".

To run the timing test of the algorithms, line 25 in [main.cpp](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_1/main.cpp) must be set to "true" like displayed here:

    // Measures duration of general and special algorithm and prints outcome if true
    bool timing = true;

If you do not want to run the timing tests, simply change "true" to "false".

To link and compile the code, we use a [makefile](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_1/makefile). In a Linux command line, write the following:

    make all

This will produce a file called main.exe. To run it, write the following in a Linux command line:

    ./main.exe n
    
where n is the number of grid points. 




