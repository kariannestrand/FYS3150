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

    
    

    

