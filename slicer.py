import numpy as np
from scipy.optimize import lsq_linear


def angleWrapTo(x, y):
    '''Wraps the input array x by adding an integer multiple of 2*pi to each
    element in it such that x-y <= pi.'''
    return x - 2*np.pi*np.rint((x-y)/(2.*np.pi))


def getPhi(nSlices):
    '''Outputs the longitudes of the center of each slice in radians for a
    given number of slices.'''
    dphi = 2*np.pi/nSlices
    return np.linspace(-np.pi, np.pi, nSlices+1)[:-1] + dphi/2


def getG(phiObs, nSlices, phiIn=None):
    '''Computes the matrix G, which relates the slice albedo
    to the phase curve.  If phiIn is None, reverts to measuring raw brightness
    without respect to an illumination source.

    args:
        phiObs: An array of floats containing the latitude over which each
            measurement was taken.
        nSlices: The number of slices for the model to use.
        phiIn (optional): An array of floats containing the latitude over which
            the illumination source (the parent star) was.
    '''
    # Phi is the location of the center of each slice
    phi = getPhi(nSlices)[None, :]
    # Get the bounds of integration
    phiPlus = phi + np.pi/nSlices
    phiMinus = phi - np.pi/nSlices
    # Clip bounds to within the visible area
    phiObs = angleWrapTo(phiObs[:, None], phi)
    phiPlus = np.clip(phiPlus, phiObs-np.pi/2, phiObs+np.pi/2)
    phiMinus = np.clip(phiMinus, phiObs-np.pi/2, phiObs+np.pi/2)
    if phiIn is not None:
        # Clip bounds also to within illuminated area
        phiIn = angleWrapTo(phiIn[:, None], phi)
        phiPlus = np.clip(phiPlus, phiIn-np.pi/2, phiIn+np.pi/2)
        phiMinus = np.clip(phiMinus, phiIn-np.pi/2, phiIn+np.pi/2)
        # Evaluate the upper and lower bounds of the integral
        return (.5*phiPlus*np.cos(phiIn-phiObs) - .25*np.sin(phiIn+phiObs-2*phiPlus) +
                -.5*phiMinus*np.cos(phiIn-phiObs) + .25*np.sin(phiIn+phiObs-2*phiMinus))
    return np.sin(phiPlus - phiObs) - np.sin(phiMinus - phiObs)


def toPhaseCurve(phiObs, j, phiIn=None, jErr=None, residNoise=0., G=None):
    '''Computes the phase curve from slice brightnesses, with optional uncertainties'''
    # Compute G if it was not provided
    if G is None:
        G = getG(phiObs, len(j), phiIn)
    mean = np.matmul(G, j)
    if jErr is not None:
        # If jErr is standard deviations, get the covariance matrix equivalent
        if jErr.ndim == 1:
            jErr = np.diag(jErr**2)
        # Make predictions for the requested points
        err = np.matmul(G, np.matmul(jErr + residNoise**2 * np.eye(len(jErr)), G.T))
        return mean, err
    return mean


def fromPhaseCurve(phiObs, flux, nSlices, phiIn=None, priorStd=1e3, G=None, brightnessMin=0,
                   brightnessMax=np.inf, fullOutput=False):
    '''Computes the slice brightnesses from a given phase curve and prior uncertainty.'''
    # Compute G if it was not provided
    if G is None:
        G = getG(phiObs, nSlices, phiIn)
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
