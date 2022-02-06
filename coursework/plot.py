from matplotlib import pyplot as plt
import numpy as np


def plot_results(A, x, b):
    n = len(A)

    res = np.dot(A, x)

    res_inf = np.array([elem.lower for elem in res])
    res_sup = np.array([elem.upper for elem in res])

    b_inf = np.array([elem.lower for elem in b])
    b_sup = np.array([elem.upper for elem in b])

    plt.figure(figsize=(18, 8))
    plt.plot(b_inf, label='b_inf')
    plt.plot(res_inf, label='Ax_inf', ls='--')
    plt.plot(b_sup, label='b_sup')
    plt.plot(res_sup, label='Ax_sup', ls='--')
    plt.title('Results')
    plt.ylabel('value')
    plt.grid()
    plt.legend()
    plt.show()
