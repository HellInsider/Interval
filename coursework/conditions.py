import numpy as np
from kaucherpy import *


def check_osob_matrix(A_inf, A_sup):
    n = len(A_inf)
    A_mid = np.array([[(A_inf[i, j] + A_sup[i, j]) / 2 for j in range(n)] for i in range(n)])

    if np.linalg.det(A_mid) == 0:
        return True

    A_rad = np.array([[(-A_inf[i, j] + A_sup[i, j]) for j in range(n)] for i in range(n)])
    diag = np.array([np.dot(A_rad, np.fabs(np.linalg.inv(A_mid)))[i, i] for i in range(n)])

    if np.max(diag) >= 1:
        return True

    return False


def check_sub_grad(A_inf, A_sup):
    n = len(A_inf)

    for i in range(n):
        for j in range(n):
            if A_inf[i][j] > A_sup[i][j]:
                num = A_inf[i][j]
                A_inf[i][j] = A_sup[i][j]
                A_sup[i][j] = num

    C_plus = [[Kaucher(0, 0) for j in range(n)] for i in range(n)]
    C_minus = [[Kaucher(0, 0) for j in range(n)] for i in range(n)]

    for i in range(n):
        for j in range(n):
            C_plus[i][j] = Kaucher(max(0, A_inf[i][j]), max(0, A_sup[i][j]))
            C_minus[i][j] = Kaucher(min(0, A_inf[i][j]), min(0, A_sup[i][j]))

    opC = [[Kaucher(0, 0) for j in range(2 * n)] for i in range(2 * n)]
    for i in range(n):
        for j in range(n):
            opC[i][j] = C_plus[i][j]
            opC[i + n][j + n] = C_plus[i][j]
            opC[i + n][j] = C_minus[i][j]
            opC[i][j + n] = C_minus[i][j]

    opc_inf = np.array([[opC[i][j].lower for j in range(n)] for i in range(n)])
    opc_sup = np.array([[opC[i][j].upper for j in range(n)] for i in range(n)])

    return not check_osob_matrix(opc_inf, opc_sup)
