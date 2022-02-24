import matplotlib.pyplot as plt
import numpy as np

A = np.loadtxt("matrix_n_phi_1.txt")
plt.xlim(0, len(A[0]))
plt.ylim(0, len(A))
plt.imshow(A, cmap='hot', interpolation='nearest', aspect='auto')
plt.show()

plt.savefig("heatmatrix.phg")