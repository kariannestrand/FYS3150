import numpy as np
import matplotlib.pyplot as plt

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

    Rx = np.array(Rx)

    with open(filename) as f:
        for line in f.read().split("\n")[1::4]:
            cols = line.split()
            if len(cols) > 1:
                Ry.append(cols[i])

    Ry = np.array(Ry)

    with open(filename) as f:
        for line in f.read().split("\n")[2::4]:
            cols = line.split()
            if len(cols) > 1:
                Rz.append(cols[i])

    Rz = np.array(Rz)
    t = np.linspace(0,100,len(Rz))
