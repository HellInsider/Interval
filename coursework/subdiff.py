import numpy as np
from kaucherpy import *


def sti_dot_mtx(matrix):
    pos = matrix.copy()
    neg = matrix.copy()
    pos[pos < 0] = 0
    neg[neg > 0] = 0
    neg = np.fabs(neg)
    return np.block([[pos, neg], [neg, pos]])


def back_sti_vector(sti_vec):
    mid_idx = sti_vec.shape[0] // 2
    return sti_vec[:mid_idx], sti_vec[mid_idx:]


def sti_vector(vec_inf, vec_sup):
    return np.append(vec_inf, vec_sup)


def compute_sub_grad2(D, i, j, g, h):
    n = D.shape[0] // 2;
    ainf = g.lower
    asup = g.upper
    b_inf = h.lower
    b_sup = h.upper

    k = 0
    m = 0
    if ainf * asup > 0:
        k = 0 if ainf > 0 else 2
    else:
        k = 1 if ainf < asup else 3

    if b_inf * b_sup > 0:
        m = 1 if b_inf > 0 else 3
    else:
        m = 2 if b_inf <= b_sup else 4

    case = 4 * k + m
    if case == 1:
        D[i, j] = ainf
        D[i + n, j + n] = asup
    elif case == 2:
        D[i, j] = asup
        D[i + n, j + n] = asup
    elif case == 3:
        D[i, j] = asup
        D[i + n, j + n] = ainf
    elif case == 4:
        D[i, j] = ainf
        D[i + n, j + n] = ainf
    elif case == 5:
        D[i, j + n] = ainf
        D[i + n, j + n] = asup
    elif case == 6:
        if ainf * b_sup < asup * b_inf:
            D[i, j + n] = ainf
        else:
            D[i, j] = asup
        if ainf * b_inf > asup * b_sup:
            D[i + n, j] = ainf
        else:
            D[i + n, j + n] = asup
    elif case == 7:
        D[i, j] = asup
        D[i + n, j] = ainf
    elif case == 8:
        pass
    elif case == 9:
        D[i, j + n] = ainf
        D[i + n, j] = asup
    elif case == 10:
        D[i, j + n] = ainf
        D[i + n, j] = ainf
    elif case == 11:
        D[i, j + n] = asup
        D[i + n, j] = ainf
    elif case == 12:
        D[i, j + n] = asup
        D[i + n, j] = asup
    elif case == 13:
        D[i, j] = ainf
        D[i + n, j] = asup
    elif case == 14:
        pass
    elif case == 15:
        D[i, j + n] = asup
        D[i + n, j + n] = ainf
    elif case == 16:
        if ainf * b_inf > asup * b_sup:
            D[i, j] = ainf
        else:
            D[i, j + n] = -asup
        if ainf * b_sup < asup * b_inf:
            D[i + n, j + n] = ainf
        else:
            D[i + n, j] = asup
    return D


def sub_diff_method2(A, b, lr=1, acc=1e-6, max_iter=100):
    n = len(A)
    A_mid = np.array([[(A[i][j].lower + A[i][j].upper) / 2 for j in range(n)] for i in range(n)])
    b_inf = np.array([elem.lower for elem in b])
    b_sup = np.array([elem.upper for elem in b])
    sti_vec = sti_vector(b_inf, b_sup)
    A_block = sti_dot_mtx(A_mid)
    cur_x = np.dot(np.linalg.inv(A_block), sti_vec)
    prev_x = cur_x
    started = False
    cur_iter = 0
    worklist = [cur_x]
    while not started or np.linalg.norm(cur_x - prev_x) > acc:
        started = True
        if cur_iter > max_iter:
            # print("Too many iterations")
            cur_iter -= 1
            break
        prev_x = [elem for elem in cur_x]
        F = np.zeros((2 * n, 2 * n))
        for i in range(n):
            s = Kaucher(0, 0)
            for j in range(n):
                g = A[i][j]
                h = Kaucher(prev_x[j], prev_x[j + n])
                t = g * h
                s = s + t
                F = compute_sub_grad2(F, i, j, g, h)
            t = int_sub(s, b[i])
            cur_x[i] = t.lower
            cur_x[i + n] = t.upper
        dx = np.linalg.solve(F, cur_x)
        cur_x = prev_x - lr * dx
        worklist.append(cur_x)
        cur_iter += 1
        x = np.array(back_sti_vector(cur_x))
        return x.T, cur_iter


def compute_sum1(H, x, i):
    return sum([H[i][j] * x[j] for j in range(i)])


def compute_sum2(G, x, i):
    return sum([G[i][j] * x[j] for j in range(i + 1, len(G))])


def back_mul(op1, op2):
    return op1 * Kaucher(1 / op2.lower, 1 / op2.upper)


def int_sub(op1, op2):
    return Kaucher(op1.lower - op2.lower, op1.upper - op2.upper)
