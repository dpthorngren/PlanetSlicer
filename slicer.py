import numpy as np


def angleWrap(x):
    '''Wraps input angles x to [-pi,pi)'''
    return (x+np.pi) % (2*np.pi) - np.pi


def getPhi(nSlices):
    dphi = 2*np.pi/nSlices
    return np.linspace(-np.pi, np.pi, nSlices+1)[:-1] + dphi/2


def getG(xi, nSlices):
    '''Computes the matrix G, which relates the slice brighntess
    to the phase curve.'''
    dphi = 2*np.pi/nSlices
    phi = getPhi(nSlices)
    alpha = angleWrap(xi[:, None]+phi[None, :])
    alpha_plus = np.clip(alpha+dphi/2, -np.pi/2, np.pi/2)
    alpha_minus = np.clip(alpha-dphi/2, -np.pi/2, np.pi/2)
    return np.sin(alpha_plus) - np.sin(alpha_minus)


def toPhaseCurve(j, xi, G=None):
    # Compute G if it was not provided
    if G is None:
        G = getG(xi, len(j))
    return np.matmul(G, j)
