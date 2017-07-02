import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt


def run_dpc():
    # Setup
    z = np.random.normal(0, 1, (10, 100))
    k1 = 5
    k2 = 3
    k = k1 + k2
    m, T = z.shape
    f_init = np.random.normal(0, 1, (T + k, 1))
    alpha_init = np.nan
    beta_init = np.nan

    # Train
    f = f_init
    alpha = alpha_init
    beta = beta_init
    train_iterations = 100
    loss_values = np.zeros(train_iterations, 1)

    for train_iter in range(train_iterations):
        f, alpha, beta = run_train_step(z, k, f, alpha, beta)
        loss_values[train_iter] = loss(z, k, f, alpha, beta)

    # Plot
    plt.plot(loss_values)
    plt.plot(f)


def run_train_step(z, k, f, alpha, beta):
    m, T = z.shape
    alpha_new, beta_new = alpha_beta(z, f, k)
    f_star = f_alpha_beta(z, k, alpha_new, beta_new)
    f_centered = f_star - np.mean(f_star)
    f_new = np.sqrt(T _ k) * f_centered / la.norm(f_centered, 2)
    return f_new, alpha_new, beta_new


def loss(z, k, f, alpha, beta):
    return mean_squared_error(z, k, f, alpha, beta)


def C_alpha(z, alpha, k):
    m, T = z.shape
    C = np.zero((m, T + k, k + 1))
    for j in range(m):
        for t in range(C.shape[1]):
            for q in range(C.shape[2]):
                if (q >= max(1, t - T + 1)) and (q <= min(k + 1, t)):
                    C[j, t, q] = z[j, t - q + 1] - alpha[j]
    return C


def f_alpha_beta(z, k, alpha, beta):
    m, T = z.shape
    D_3d = D_beta(z, beta, k)
    D = np.squeeze(np.sum(D_3d, 1))
    C_3d = C_alpha(z, alpha, k)
    sum_Cj_betaj = 0
    for j in range(m):
        sum_Cj_betaj = sum_Cj_betaj + np.squeeze(C_3d[j, :, :])*beta(j, :).T
    return la.solve(D, sum_Cj_betaj)


def F_f(z, f, k):
    m, T = z.shape
    F = np.zeros(T, k + 2)

    for t in range(F.shape[0]):
        F[t, :] = [f(t: (t + k)).T 1]


def alpha_beta(z, f, k):
    m, T = z.shape
    F = F_f(z, f, k)
    FtF_inv - la.inv(F.T * F)
    FtF_inv_Ft = FtF_inv * F.T
    tmp = z * FtF_inv_Ft.T
    alpha = tmp(:, k+2).T
    beta = tmp(:, : (k+1))


def mean_squared_error(z, k, f, alpha, beta):
    m, T = z.shape
    sum_squared_error = 0
    for j in range(m):
        for t in range(T):
            z_jt_predict = alpha[j]
            for i in range(k):
                z_jt_predict = z_jt_predict + beta[j, i + 1] * f(t + i)
            sum_squared_error = sum_squared_error + (z[j, t] - z_jt_predict)**2
    return sum_squared_error / (T * m)
