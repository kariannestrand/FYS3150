import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d



#one particle

filename = "Euler_r.txt"
Rx = []
Ry = []
Rz = []
with open(filename) as f:
        for line in f.read().split("\n")[0::4]:
            Rx.append(line)


Rx = np.array([float(x) for x in Rx[:-1]])

with open(filename) as f:
        for line in f.read().split("\n")[1::4]:
            Ry.append(line)

Ry = np.array([float(x) for x in Ry])

with open(filename) as f:
        for line in f.read().split("\n")[2::4]:
            Rz.append(line)

Rz = np.array([float(x) for x in Rz])


fig = plt.figure()
ax = plt.axes(projection="3d")
ax.plot3D(Rx, Ry, Rz, 'gray')

plt.savefig('3d.pdf')
plt.show()


t = np.linspace(0, 100, len(Rz))

q = 1.
V0 = 9.65*10e8
m = 40.078
d = 10e4
omega_z = np.sqrt(2*q*V0/(m*d**2))
print(omega_z)


fig = plt.figure()
plt.plot(Rx, Ry)
plt.savefig('xy_movement.pdf')
plt.show()

fig = plt.figure()
plt.plot(t, Rz)
plt.savefig('z_movement.pdf')
plt.show()




#several particles
n_particles = 2 # number of particles

for i in range(n_particles):
    filename = "Euler_r.txt"
    Rx = []
    Ry = []
    Rz = []
    with open(filename) as f:
        for line in f.read().split("\n")[0::4]:
            cols = line.split()
            if len(cols) > 1:
                Rx.append(cols[i])
            else:
                Rx.append(cols[0])

    Rx = np.array(Rx)

    with open(filename) as f:
        for line in f.read().split("\n")[1::4]:
            cols = line.split()
            if len(cols) > 1:
                Ry.append(cols[i])
            else:
                Ry.append(cols[0])

    Ry = np.array(Ry)

    with open(filename) as f:
        for line in f.read().split("\n")[2::4]:
            cols = line.split()
            if len(cols) > 1:
                Rz.append(cols[i])
            else:
                Rz.append(cols[0])

    Rz = np.array(Rz)




