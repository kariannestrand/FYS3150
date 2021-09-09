import numpy as np
import matplotlib.pyplot as plt

x_exact = []; u_exact = []
with open("problem2.txt") as f:
    f.readline()
    for line in f.readlines():
        line = line.split()
        x_exact.append(float(line[0])); u_exact.append(float(line[1]))

x10 = []; u10 = []
with open("problem7.txt") as f:
    f.readline()
    for line in f.readlines():
        line = line.split()
        x10.append(float(line[0])); u10.append(float(line[1]))

x100 = []; u100 = []
with open("problem7.txt") as f:
    f.readline()
    for line in f.readlines():
        line = line.split()
        x100.append(float(line[0])); u100.append(float(line[1]))

x1000 = []; u1000 = []
with open("problem7.txt") as f:
    f.readline()
    for line in f.readlines():
        line = line.split()
        x1000.append(float(line[0])); u1000.append(float(line[1]))


plt.plot(x_exact, u_exact, linestyle="-", marker="o")
plt.plot(x10, u10, inestyle="-", marker="o", color='red')
plt.plot(x100, u100, inestyle="-", marker="o", color='green')
plt.plot(x100, u100, inestyle="-", marker="o", color='yellow')
plt.title('Plot of $u(x)=1-(1-e^{-10})x-e^{-10x}$')
plt.xlabel('x-values')
plt.ylabel('u(x)')
plt.savefig("problem7.pdf")
plt.show()