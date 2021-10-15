import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

#two particles, velocity

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
ax.plot3D(Rx1, Ry1, Rz1, 'blue')
ax.plot3D(Rx2, Ry2, Rz2, 'red')

plt.savefig('3d_2particles.pdf')
plt.show()

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


'''
#several particles in loop
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

    
    

    

