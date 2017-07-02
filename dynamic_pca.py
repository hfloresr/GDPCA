import numpy as np
from numpy import linalg as la


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


def alpha_beta(z, f, k):
    m, T = z.shape
