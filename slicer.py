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
        return (2*phiPlus*np.cos(phiIn-phiObs) - np.sin(phiIn+phiObs-2*phiPlus) +
                -2*phiMinus*np.cos(phiIn-phiObs) + np.sin(phiIn+phiObs-2*phiMinus)) / 3.
    return np.sin(phiPlus - phiObs) - np.sin(phiMinus - phiObs)


def toPhaseCurve(phiObs, j, phiIn=None, jErr=None, residNoise=0., G=None):
    '''Computes the phase curve from slice brightnesses, with optional uncertainties.

    args:
        phiObs: The longitude above which each observation is to be simulated, as an array of length
            [observations], in radians from -π to π.
        j: The brightness (energy per area) for each slice.  Alternatively if `phiIn` is given,
            it is the surface albedo for each slice.  Either way, it should be an array of
            length [slices].
        phiIn (optional): The longitude above which is the source of illumination, in radians from
            -π to π.  If not specified, the object is assumed to be self-luminous.  Otherwise,
            it should be an array of shape [observations].  The default is None.
        jErr (optional): The uncertainty on j, from which to compute the uncertainty in the
            predicted phase curve.  If unspecified, then no uncertainty is computed.  Otherwise,
            it should be an array of shape [slices].  The default is None.
        residNoise (optional): Additional white noise to add to the unceretainty in the predicted
            phase curve, as a float.  Ignored if jErr is None.  The default is 0.
        G (optional): The matrix which relates slice brightness to observed total brightnesses,
            whose shape is [observations X slices].  This is calculated automatically if not
            specified.  The default is None.

    returns:
        The predicted observed flux relative to the incident flux as an array of length
        [observations].  If jErr was specified, then the uncertainty is also returned
        with the same shape.'''
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
    '''Computes the slice brightnesses from a given phase curve with optional uncertainties.

    args:
        phiObs: The longitude above which each observation was made, as an array of length
            [observations], in radians from -π to π.
        flux: The brightness observed for each observation, as an array of length [observations].
        nSlices: The number of slices to use in the model, which must be an integer greater than 0.
        phiIn (optional): The longitude above which is the source of illumination, in radians from
            -π to π.  If not specified, the object is assumed to be self-luminous.  Otherwise, it
            should be an array of shape [observations].  The default is None.
        priorStd (optional): The standard deviation of the prior on the slice brightnesses.
            Setting a large value can mimic a flat prior, but can cause numerical problems if G
            is not full rank.  The default is 1000.
        G (optional): The matrix which relates slice brightness to observed total brightnesses,
            whose shape is [observations X slices].  This is calculated automatically if not
            specified.  The default is None.
        brightnessMin (optional): The minimum allowed brightness for each slice.  The default is 0.
        brightnessMax (optional): The maximum allowed brightness for each slice.  If `phiIn` is set,
            this is the max albedo which you may wish to set to 1 (the maximum possible surface
            albedo).  The default is positive infinity.
        fullOutput (optional): If True, returns a dictionary containing uncertainty information in
            addition to the usual output.  The default is False.

    returns:
        The brightnesses inferred for the slices, as an array of length [slices].  If `phiIn` was
        set, these are surface albedos.  If `fullOutput` was set to true, the function instead
        returns a dictionary containing the brightnesses (or albedos), the posterior covariance,
        the log likelihood of the model, and the residual standard deviation.'''
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
