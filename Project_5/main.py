import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa


t = 3
for i in range(t):
    U_in = pa.cx_mat()
    U_in.load("U_in_" + str(i) + "t.bin")
print(U_in)
