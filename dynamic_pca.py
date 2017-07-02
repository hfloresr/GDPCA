#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from numpy import linalg as la
from scipy.io import loadmat
import matplotlib.pyplot as plt


def run_dpc():
    data_z = loadmat('z.mat')
    data_f = loadmat('f.mat')

    # Setup
    #z = np.random.normal(0, 1, (10, 100))
    z = data_z['z']
    k1 = 5
    k2 = 3
    k = k1 + k2
    m, T = z.shape
    #f_init = np.random.normal(0, 1, (T + k, 1))
    f_init = data_f['f']
    alpha_init = np.nan
    beta_init = np.nan

    # Train
    f = f_init
    alpha = alpha_init
    beta = beta_init
    train_iterations = 100
    loss_values = np.zeros((train_iterations, 1))

    for train_iter in range(train_iterations):
        f, alpha, beta = run_train_step(z, k, f, alpha, beta)
        loss_values[train_iter] = loss(z, k, f, alpha, beta)

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 8))
    ax1.plot(loss_values)
    ax1.set_title('Train loss (reconstruction MSE)')
    ax1.set_xlabel('iteration')

    ax2.plot(f)
    ax2.set_title('Factor')
    ax2.set_xlabel('time')
    plt.show()


def run_train_step(z, k, f, alpha, beta):
    m, T = z.shape
    alpha_new, beta_new = alpha_beta(z, f, k)
    f_star = f_alpha_beta(z, k, alpha_new, beta_new)
    f_centered = f_star - np.mean(f_star)
    f_new = np.sqrt(T + k) * f_centered / la.norm(f_centered, 2)
    return f_new, alpha_new, beta_new


def loss(z, k, f, alpha, beta):
    return mean_squared_error(z, k, f, alpha, beta)


def C_alpha(z, alpha, k):
    m, T = z.shape
    C = np.zeros((m, T + k, k + 1))
    for j in range(m):
        for t in range(1, C.shape[1] + 1):
            for q in range(1, C.shape[2] + 1):
                if (q >= max(1, t - T + 1)) and (q <= min(k + 1, t)):
                    C[j, t - 1, q - 1] = z[j, t - q] - alpha[j]
    return C


def D_beta(z, beta, k):
    m, T = z.shape
    D = np.zeros((m, T + k, T + k))
    for j in range(m):
        for t in range(1, D.shape[1] + 1):
            for q in range(1, D.shape[2] + 1):
                if (q >= max(t - k, 1)) and (q <= min(t + k, T + k)):
                    for v in range(max([1, t - k, q - k]), min([t, q, T]) + 1):
                        D[j, t - 1, q -
                            1] += np.dot(beta[j, q - v], beta[j, t - v])
    return D


def f_alpha_beta(z, k, alpha, beta):
    m, T = z.shape
    D_3d = D_beta(z, beta, k)
    D = np.squeeze(np.sum(D_3d, 0))
    C_3d = C_alpha(z, alpha, k)
    sum_Cj_betaj = 0
    for j in range(m):
        sum_Cj_betaj += np.dot(np.squeeze(C_3d[j, :, :]), beta[j, :].T)
    sum_Cj_betaj = sum_Cj_betaj[:, np.newaxis]
    return la.solve(D, sum_Cj_betaj)


def F_f(z, f, k):
    m, T = z.shape
    F = np.zeros((T, k + 2))

    for t in range(F.shape[0]):
        F[t, :] = np.column_stack((f[t: (t + k + 1)].T, 1))
    return F


def alpha_beta(z, f, k):
    m, T = z.shape
    F = F_f(z, f, k)
    FtF_inv = la.inv(np.dot(F.T, F))
    FtF_inv_Ft = np.dot(FtF_inv, F.T)
    tmp = np.dot(z, FtF_inv_Ft.T)
    alpha = tmp[:, k + 1].T
    beta = tmp[:, :(k + 1)]
    return alpha, beta


def mean_squared_error(z, k, f, alpha, beta):
    m, T = z.shape
    sum_squared_error = 0
    for j in range(m):
        for t in range(T):
            z_jt_predict = alpha[j]
            for i in range(k + 1):
                z_jt_predict += beta[j, i] * f[t + i]
            sum_squared_error = sum_squared_error + (z[j, t] - z_jt_predict)**2
    return sum_squared_error / (T * m)


if __name__ == '__main__':
    run_dpc()
