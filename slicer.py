import numpy as np
from scipy.optimize import differential_evolution


def angleWrap(x):
    '''Wraps input angles x to [-pi,pi)'''
    return (x+np.pi) % (2*np.pi) - np.pi


def getPhi(nSlices):
    '''Outputs the longitudes of the center of each slice in radians for a
    given number of slices.'''
    dphi = 2*np.pi/nSlices
    return np.linspace(-np.pi, np.pi, nSlices+1)[:-1] + dphi/2


def getG(xi, nSlices):
    '''Computes the matrix G, which relates the slice brightness
    to the phase curve.'''
    dphi = 2*np.pi/nSlices
    phi = getPhi(nSlices)
    alpha = angleWrap(xi[:, None]+phi[None, :])
    alpha_plus = np.clip(alpha+dphi/2, -np.pi/2, np.pi/2)
    alpha_minus = np.clip(alpha-dphi/2, -np.pi/2, np.pi/2)
    return np.sin(alpha_plus) - np.sin(alpha_minus)


def toPhaseCurve(j, xi, G=None):
    '''Computes the phase curve from slice brightnesses.'''
    # Compute G if it was not provided
    if G is None:
        G = getG(xi, len(j))
    return np.matmul(G, j)


def fromPhaseCurve(xi, flux, nSlices, fluxErr=None, G=None, fullOutput=False,
                   brightnessMax=10.):
    '''Computes the slice brightnesses from a given phase curve.  Technically,
    this is a linear model, but it is so poorly-behaved that this approach uses
    differential evolution to solve the equations.  This can *still* be
    unreliable for larger numbers of slices, so be cautious.'''
    # Compute G if it was not provided
    if G is None:
        G = getG(xi, nSlices)
    # If fluxErr is not provided, just use flat flux uncertainties.
    if fluxErr is None:
        fluxErr = np.ones(len(flux))

    # Solve for the slice brightnesses
    j = differential_evolution(
        lambda j: (((np.matmul(G, j)-flux)/fluxErr)**2).sum(),
        list(zip(np.zeros(nSlices), brightnessMax*np.ones(nSlices))))
    if fullOutput:
        return j
    return j.x
