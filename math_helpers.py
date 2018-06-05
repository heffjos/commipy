import numpy as np

def pearson_2d(X, Y):
    """ Pearson product-moment correlation of vectors `X` and array `Y`
    Parameters
    ----------
    X : array shape (P, N)
        One-dimensional array to correlate with every column of `Y`
    Y : array shape (Q, N)
        2D array where we correlate each column of `Y` with `X`.
    Returns
    -------
    r_xy : array shape (P,Q)
        Pearson product-moment correlation of vectors `X` and the columns of
        `Y`, with one correlation value for every column of `Y`.
    """
    mc_X = X - np.mean(X, 1, keepdims = True)
    mc_Y = Y - np.mean(Y, 1, keepdims = True)

    ssX = (mc_X**2).sum(1)
    ssY = (mc_Y**2).sum(1)

    return np.dot(mc_X, mc_Y.T) / np.sqrt(np.dot(ssX[:, np.newaxis], ssY[np.newaxis]))
