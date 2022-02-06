from conditions import check_sub_grad, check_osob_matrix
from subdiff import sub_diff_method2
import numpy as np
from kaucherpy import *

###          TEST1              ###
A_inf = np.array([[2, -2],
                  [-1, 2]])
A_sup = np.array([[4, 1],
                  [2, 4]])

b_inf = np.array([-2, -2])
b_sup = np.array([2, 2])

n = len(A_inf)
a = [[Kaucher(A_inf[i, j], A_sup[i, j]) for j in range(n)] for i in range(n)]
b = [Kaucher(b_inf[i], b_sup[i]) for i in range(n)]

check_osob_matrix(A_inf, A_sup)
check_sub_grad(A_inf, A_sup)

x, iter = sub_diff_method2(a, b, lr=1)
print(f'''x = {x} \niter = {iter}''')

###          TEST2              ###
A_inf = np.array([[2, -5, -2],
                  [-3, 5, 4],
                  [-1, -2, -7]])

A_sup = np.array([[4, -1, 3],
                  [1, 7, 6],
                  [1, 1, -2]])

n = len(A_inf)

b_inf = np.array([-28, -60, -11])

b_sup = np.array([43, 29, 39])

a = [[Kaucher(A_inf[i, j], A_sup[i, j]) for j in range(n)] for i in range(n)]
b = [Kaucher(b_inf[i], b_sup[i]) for i in range(n)]

check_osob_matrix(A_inf, A_sup)
check_sub_grad(A_inf, A_sup)

x, iter = sub_diff_method2(a, b, lr=1)
print(f'''x = {x} \niter = {iter}''')

###          TEST2              ###
A_inf = np.array([[3, -5, -2],
                  [-3, 6, 5],
                  [-1, -1, -4]])

A_sup = np.array([[4, -2, 2],
                  [-1, 7, 6],
                  [0, 1, 1]])

n = len(A_inf)

b_inf = np.array([-28, -60, -11])

b_sup = np.array([43, 69, 39])

a = [[Kaucher(A_inf[i, j], A_sup[i, j]) for j in range(n)] for i in range(n)]
b = [Kaucher(b_inf[i], b_sup[i]) for i in range(n)]

check_osob_matrix(A_inf, A_sup)
check_sub_grad(A_inf, A_sup)

x, iter = sub_diff_method2(a, b, lr=1)
print(f'''x = {x} \niter = {iter}''')

###          EXPERIMENT        ###
a_new = [[Kaucher(A_inf[i, j], A_sup[i, j]) for j in range(n)] for i in range(n)]

for num in np.arange(0, 10, 0.5):
    a_new[2][2] = Kaucher(num, A_sup[2, 2])
    x, iter = sub_diff_method2(a_new, b, lr=0.5)
    if iter < 100:
        print(f'Change a33 inf to {num} \t => converged')
    else:
        print(f'Change a33 inf to {num} \t => diverged')
