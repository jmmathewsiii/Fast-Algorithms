import numpy as np


def generate_random_matrix(N, seed=None):

    rng = np.random.default_rng(seed)
    mat = rng.standard_normal((N, N))

    return mat


def sample_bivariate_function(f, x_min, x_max, y_min, y_max, N):

    x_vals = np.linspace(x_min, x_max, N)
    y_vals = np.linspace(y_min, y_max, N)

    x_mat, y_mat = np.meshgrid(x_vals, y_vals, indexing='ij')

    result = f(x_mat, y_mat)

    return result


def generate_singular_values(N):


