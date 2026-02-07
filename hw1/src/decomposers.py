import numpy as np


def get_qr_decomposition(matrix):

    Q, R = np.linalg.qr(matrix)

    return Q, R


def get_svd(matrix):

    U, S, Vt = np.linalg.svd(matrix)

    return U, S, Vt
