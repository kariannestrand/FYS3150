import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

n = np.linspace(10,100,10)
t = np.array([161, 729, 1715, 3079, 4861, 6995, 9530, 12498, 15746, 19624])

reg = linregress(np.log10(n),np.log10(t))
print(reg)

plt.title('Scaling of number of required transformations with matrix size N')
plt.plot(np.log10(n), np.log10(t), label = 'Slope = 2.07')
plt.legend()
plt.xlabel('log10(N)', size = 10)
plt.ylabel('log10(T)', size = 10)
plt.tight_layout()
plt.savefig('scaling.pdf')
plt.show()