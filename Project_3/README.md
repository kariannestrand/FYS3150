## Project 3 - The Penning Trap

The goal of this project is to numerically simulate a Penning trap, and exploring different aspects of it. To solve the ordinary differential equations, we have implemented the Runge-Kutta 4, and the Euler Cromer method. We made a Particle class to create the particles, and a Penning Trap class to simulate the trap itself. The classes both have a header file (.hpp) to define all class methods and variables needed, and a source file (.cpp) to implement all methods used.

In main.cpp we have set all the initial conditions needed to simulate the Penning trap. There are also some different aspects of the trap that is possible to turn on and off, by changing between the options True and False. These aspects are particle interactions, modification by applying a timevarying potential, and also to choose between  running the Euler Cromer or the RK4 method. 

To build and run the code, we have made a makefile that builds and runs the code by writing MAKE in a Linux or Unix command line. 
