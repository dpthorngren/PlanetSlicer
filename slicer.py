import numpy as np
from scipy.optimize import lsq_linear


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


def toPhaseCurve(xi, j, jErr=None, residNoise=0., G=None):
    '''Computes the phase curve from slice brightnesses, with optional uncertainties'''
    # Compute G if it was not provided
    if G is None:
        G = getG(xi, len(j))
    mean = np.matmul(G, j)
    if jErr is not None:
        # If jErr is standard deviations, get the covariance matrix equivalent
        if jErr.ndim == 1:
            jErr = np.diag(jErr**2)
        # Make predictions for the requested points
        err = np.matmul(G, np.matmul(jErr + residNoise, G.T))
        return mean, err
    return mean


def fromPhaseCurve(xi, flux, nSlices, priorStd=1e3, G=None, brightnessMin=0,
                   brightnessMax=np.inf, fullOutput=False):
    '''Computes the slice brightnesses from a given phase curve and prior uncertainty.'''
    # Compute G if it was not provided
    if G is None:
        G = getG(xi, nSlices)
    # Compute with bounded ridge regression using the pseudo-observations method.
    X = np.vstack([G, np.eye(nSlices)/priorStd])
    y = np.concatenate([flux, np.mean(flux) * np.ones(nSlices) / (2*priorStd)])
    f = lsq_linear(X, y, bounds=[brightnessMin, brightnessMax])
    if not fullOutput:
        return f.x
    # Compute the log posterior probability
    sigma = np.sqrt(2.*f.cost/(len(flux)-nSlices))
    logLike = -len(flux)*np.log(2*np.pi)/2. - len(flux)*sigma - f.cost / sigma**2
    # Get the uncertainties on the brightness
    errors = sigma**2 * np.linalg.inv(np.matmul(X.T, X))
    return {'brightness': f.x, 'brightnessCov': errors, 'logLike': logLike, 'residStd': sigma}
