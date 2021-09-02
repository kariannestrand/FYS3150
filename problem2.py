import numpy as np
import matplotlib.pyplot as plt

x = []; u = []
with open("problem2.txt") as f:
    f.readline()
    for line in f.readlines():
        line = line.split()
        x.append(float(line[0])); u.append(float(line[1]))


plt.plot(x,u, linestyle="", marker="o")
plt.title('Plot of $u(x)=1-(1-e^{-10})x-e^{-10x}$')
plt.xlabel('x-values')
plt.ylabel('u(x)')
#plt.savefig("problem2.pdf")
plt.show()