import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


#two particles, velocity

filename = "RK4_v.txt"
Vx1 = []
Vy1 = []
Vz1 = []
with open(filename) as f:
        for line in f.read().split("\n")[0::4]:
            particles = line.split()
            if particles:
                Vx1.append(particles[0])

Vx1 = np.array([float(x) for x in Vx1])

with open(filename) as f:
        for line in f.read().split("\n")[1::4]:
            particles = line.split()
            if particles:
                Vy1.append(particles[0])

Vy1 = np.array([float(x) for x in Vy1])

with open(filename) as f:
        for line in f.read().split("\n")[2::4]:
            particles = line.split()
            if particles:
                Vz1.append(particles[0])

Vz1 = np.array([float(x) for x in Vz1])

Vx2 = []
Vy2 = []
Vz2 = []
with open(filename) as f:
        for line in f.read().split("\n")[0::4]:
            particles = line.split()
            if particles:
                Vx2.append(particles[1])


Vx2 = np.array([float(x) for x in Vx2])

with open(filename) as f:
        for line in f.read().split("\n")[1::4]:
            particles = line.split()
            if particles:
                Vy2.append(particles[1])

Vy2 = np.array([float(x) for x in Vy2])

with open(filename) as f:
        for line in f.read().split("\n")[2::4]:
            particles = line.split()
            if particles:
                Vz2.append(particles[1])

Vz2 = np.array([float(x) for x in Vz2])

#two particles, position

filename = "RK4_r.txt"
Rx1 = []
Ry1 = []
Rz1 = []
with open(filename) as f:
        for line in f.read().split("\n")[0::4]:
            particles = line.split()
            if particles:
                Rx1.append(particles[0])

Rx1 = np.array([float(x) for x in Rx1])

with open(filename) as f:
        for line in f.read().split("\n")[1::4]:
            particles = line.split()
            if particles:
                Ry1.append(particles[0])

Ry1 = np.array([float(x) for x in Ry1])

with open(filename) as f:
        for line in f.read().split("\n")[2::4]:
            particles = line.split()
            if particles:
                Rz1.append(particles[0])

Rz1 = np.array([float(x) for x in Rz1])

Rx2 = []
Ry2 = []
Rz2 = []
with open(filename) as f:
        for line in f.read().split("\n")[0::4]:
            particles = line.split()
            if particles:
                Rx2.append(particles[1])


Rx2 = np.array([float(x) for x in Rx2])

with open(filename) as f:
        for line in f.read().split("\n")[1::4]:
            particles = line.split()
            if particles:
                Ry2.append(particles[1])

Ry2 = np.array([float(x) for x in Ry2])

with open(filename) as f:
        for line in f.read().split("\n")[2::4]:
            particles = line.split()
            if particles:
                Rz2.append(particles[1])

Rz2 = np.array([float(x) for x in Rz2])


# plotting for 1 and 2 particles in problem 9

fig = plt.figure()
ax = plt.axes(projection="3d")
plt.tight_layout()
ax.plot3D(Rx1, Ry1, Rz1, 'blue', label='Trajectory of particle 1 with interactions')
ax.plot3D(Rx2, Ry2, Rz2, 'red', label ='Trajectory of particle 2 with interactions')
plt.legend()

plt.savefig('3d_2particles.pdf')
plt.show()

'''
fig = plt.figure()
plt.plot(Rx1, Ry1)
plt.plot(Rx2, Ry2)
plt.savefig('xy_movement_2part.pdf')
plt.show()

t = np.linspace(0,100,len(Rz1))
fig = plt.figure()
plt.plot(t, Rz1)
plt.savefig('z_movement.pdf')
plt.show()

fig = plt.figure()
plt.plot(Rx1, Vx1)
plt.plot(Rx2, Vx2)
plt.savefig('phase_space_x.pdf')
plt.show()

fig = plt.figure()
plt.plot(Ry1, Vy1)
plt.plot(Ry2, Vy2)
plt.savefig('phase_space_y.pdf')
plt.show()

fig = plt.figure()
plt.plot(Rz1, Vz1)
plt.plot(Rz2, Vz2)
plt.savefig('phase_space_z.pdf')
plt.show()


#everything below here not working

# plotting relative errors, have to import the data from file in the right way 

relative_error_euler = False
if relative_error_euler:

    # function calculating relative error:
    def rel_error(exact, euler):
        return np.abs((exact - euler)/euler)

    n = [1, 2, 3, 4, 5]  #change this to the value of dt in the different cases
    for i in range(len(n)):
        # reading from files, extracting x, u from exact_n.txt and v from approx_n.txt:
        filename_exact = "exact_" + str(n[i]) + ".txt"
        #import the 5 different exact solutions here
        
        filename_approx = "euler_" + str(n[i]) + ".txt"
        #import the 5 different euler solutions here

        filename_approx = "time_" + str(n[i]) + ".txt"
        #import the 5 different time-arrays here

        # plotting relative error by calling function rel_error with euler and exact:
        plt.plot(t, rel_error(exact, euler), "--", label = "dt = " + str(n[i]))
        plt.title("Relative error", fontsize = 20)
        plt.xlabel("t", size = 20); plt.xticks(size = 15)
        plt.ylabel("$\epsilon$", size = 20); plt.yticks(size = 15)
        plt.legend(loc = "upper left", fontsize = 10)
        plt.tight_layout()
        plt.savefig("relative_error_euler.pdf")
    plt.show()


relative_error_rk4 = False
if relative_error_rk4:

    # function calculating relative error:
    def rel_error(exact, rk4):
        return np.abs((exact - rk4)/rk4)

    n = [1, 2, 3, 4, 5]  #change this to the value of dt in the different cases
    for i in range(len(n)):
        # reading from files, extracting x, u from exact_n.txt and v from approx_n.txt:
        filename_exact = "exact_" + str(n[i]) + ".txt"
        #import the 5 different exact solutions here
        
        filename_approx = "rk4_" + str(n[i]) + ".txt"
        #import the 5 different rk4 solutions here

        filename_approx = "time_" + str(n[i]) + ".txt"
        #import the 5 different time-arrays here  

        # plotting relative error by calling function rel_error with exact and rk4:
        plt.plot(t, rel_error(exact, rk4), "--", label = "dt = " + str(n[i]))
        plt.title("Relative error", fontsize = 20)
        plt.xlabel("t", size = 20); plt.xticks(size = 15)
        plt.ylabel("$\epsilon$", size = 20); plt.yticks(size = 15)
        plt.legend(loc = "upper left", fontsize = 10)
        plt.tight_layout()
        plt.savefig("relative_error_rk4.pdf")
    plt.show()

#several particles in loop, not working
n_particles = 2 # number of particles

for i in range(n_particles):
    filename = "RK4_r.txt"
    Rx = []
    Ry = []
    Rz = []
    with open(filename) as f:
        for line in f.read().split("\n")[0::4]:
            particles = line.split() 
            if len(particles) > 1:   
                Rx.append(particles[i])
                
    
    Rx = np.array([float(x) for x in Rx])
    
    with open(filename) as f:
        for line in f.read().split("\n")[0::4]:
            particles = line.split() 
            if len(particles) > 1:   
                Ry.append(particles[i])

    Ry = np.array([float(x) for x in Ry])

    with open(filename) as f:
        for line in f.read().split("\n")[0::4]:
            particles = line.split() 
            if len(particles) > 1:   
                Rz.append(particles[i])

    Rz = np.array([float(x) for x in Rz])
    print(Rz)

'''
    

    

