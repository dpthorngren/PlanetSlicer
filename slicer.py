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


def toPhaseCurve(j, xi, G=None):
    '''Computes the phase curve from slice brightnesses.'''
    # Compute G if it was not provided
    if G is None:
        G = getG(xi, len(j))
    return np.matmul(G, j)


def fromPhaseCurve(xi, flux, nSlices, priorStd=1e3, xiPredict=None, G=None,
                   brightnessMin=0, brightnessMax=np.inf, fullOutput=False):
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
    logProb = -len(flux)*np.log(2*np.pi)/2. - len(flux)*sigma - f.cost / sigma**2
    # Get the uncertainties on the brightness
    errors = sigma**2 * np.linalg.inv(np.matmul(X.T, X))
    output = {'brightness': f.x, 'brightnessCov': errors, 'logProb': logProb}
    if xiPredict is not None:
        # Make predictions for the requested points
        gPredict = getG(xi, nSlices)
        output['fluxPredict'] = np.matmul(gPredict, f.x)
        output['fluxPredictCov'] = \
            np.matmul(gPredict, np.matmul(errors, gPredict.T)) + sigma**2 * np.eye(len(xi))
    return output
