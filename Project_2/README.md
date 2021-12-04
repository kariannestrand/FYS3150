# Eigenvalue problems

We have organized the code for the project in two source files and one header file, [jacobi_rotation.hpp](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_2/jacobi_rotation.hpp), [jacobi_rotation.cpp](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_2/jacobi_rotation.cpp) and 
[main.cpp](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_2/main.cpp). 
[jacobi_rotation.hpp](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_2/jacobi_rotation.hpp) declares all the methods used for the project, in the class called jacobi_rotation. 
[jacobi_rotation.cpp](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_2/jacobi_rotation.cpp) contains all the calculations. [main.cpp](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_2/main.cpp) defines all the initial values, writes the necessary output to file and calls all necessary functions/methods.

To check that the numerical results agree with the analytical result, make sure that line 56 in main is set to "true" to print the output, like displayed here:

    bool print_eig = true;

If you do not want to print the output, change "true" to "false".

To get the output necessary for the plots written to file, line 89 in [main.cpp](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_2/main.cpp) must be set to "true" like displayed here:

    bool make_txt_files = true;

If you do not want the output written to file, change "true" to "false".

To print out the number of required transformations, set line 80 to "true" like displayed here:

    bool print_count = true;


To link and compile the code, we use a [makefile](https://github.com/mariaoftedahl/FYS3150/blob/main/Project_2/makefile). In a command line in Linux or Unix, write the following:

    make all

This will produce a file called main.exe. To run it, write the following in a command line:

    ./main.exe n
    
where n=N-1 where N is the size of the matrix.
