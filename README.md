# PlanetSlicer
This is a small package for fitting brightness maps to phase curves using the "orange-slice" method.  For self-luminous objects, this follows [Cowan and Agol (2008) "Inverting Phase Functions to Map Exoplanets"](https://doi.org/10.1017/S1743921308027166).  We have extended this approach to diffuse reflected light assuming [Lambertian reflectance](https://en.wikipedia.org/wiki/Lambertian_reflectance), in Mayorga et. al (in prep).  In both cases, the model supposes that a spherical object can be divided into slices of constant brightness (or albedo) which may be integrated to yield the total flux observed, given the angle(s) of observation `phiObs`.

## Installation
The package consists only of a single file, slicer.py, which may be added to the Python path or copied into the directory where it is to be used.  The only requirements are [Numpy](https://numpy.org/) and [Scipy](https://www.scipy.org/), which may be installed manually or using [pip](https://pypi.org/project/pip/).

## Usage
This package contains two key functions: toPhaseCurve and fromPhaseCurve (their docstrings are reproduced below).  The former integrates the brightness for each slice to calculate the observed total flux from the object, given the longitude of observation.  The latter does the opposite, estimating the brightness of the slices from a set of observed total flux (the phase curve).

If `phiIn` is specified for either of these functions, then instead of modeling a self-luminous object, they model one diffusely reflecting light, where `phiIn` gives the longitude pointing towards the source of the light, presumably the parent star.  In this case, the output of fromPhaseCurve and input of toPhaseCurve is no longer the slice brightness, but instead its surface albedo.  Phi (in `phiIn` and `phiOut`) always refers to a planet longitude in radians from -π to π; the zero point (prime meridian) may be arbitrarily chosen, so long as it is consistently used.


### toPhaseCurve()
Computes the phase curve from slice brightnesses, with optional uncertainties.

Arguments:
* **phiObs**: The longitude above which each observation is to be simulated, as an array of length [observations], in radians from -π to π.
* **j**: The brightness (energy per area) for each slice.  Alternatively if `phiIn` is given, it is the surface albedo for each slice.  Either way, it should be an array of length [slices].
* **phiIn** (optional): The longitude above which is the source of illumination, in radians from -π to π.  If not specified, the object is assumed to be self-luminous.  Otherwise, it should be an array of shape [observations].  The default is None.
* **jErr** (optional): The uncertainty on j, from which to compute the uncertainty in the predicted phase curve.  If unspecified, then no uncertainty is computed.  Otherwise, it should be an array of shape [slices].  The default is None.
* **residNoise** (optional): Additional white noise to add to the unceretainty in the predicted phase curve, as a float.  Ignored if jErr is None.  The default is 0.
* **G** (optional): The matrix which relates slice brightness to observed total brightnesses, whose shape is [observations X slices].  This is calculated automatically if not specified.  The default is None.

Returns: The predicted observed flux relative to the incident flux as an array of length [observations].  If jErr was specified, then the uncertainty is also returned with the same shape.

### fromPhaseCurve()
Computes the slice brightnesses from a given phase curve with optional uncertainties.

Arguments:
* **phiObs**: The longitude above which each observation was made, as an array of length [observations], in radians from -π to π.
* **flux**: The brightness observed for each observation, as an array of length [observations].
* **nSlices**: The number of slices to use in the model, which must be an integer greater than 0.
* **phiIn** (optional): The longitude above which is the source of illumination, in radians from -π to π.  If not specified, the object is assumed to be self-luminous.  Otherwise, it should be an array of shape [observations].  The default is None.
* **priorStd** (optional): The standard deviation of the prior on the slice brightnesses.  Setting a large value can mimic a flat prior, but can cause numerical problems if G is not full rank.  The default is 1000.
* **G** (optional): The matrix which relates slice brightness to observed total brightnesses, whose shape is [observations X slices].  This is calculated automatically if not specified.  The default is None.
* **brightnessMin** (optional): The minimum allowed brightness for each slice.  The default is 0.
* **brightnessMax** (optional): The maximum allowed brightness for each slice.  If `phiIn` is set, this is the max albedo which you may wish to set to 1 (the maximum possible surface albedo).  The default is positive infinity.
* **fullOutput** (optional): If True, returns a dictionary containing uncertainty information in addition to the usual output.  The default is False.

Returns: The brightnesses inferred for the slices, as an array of length [slices].  If `phiIn` was set, these are surface albedos.  If `fullOutput` was set to true, the function instead returns a dictionary containing the brightnesses (or albedos), the posterior covariance, the log likelihood of the model, and the residual standard deviation.
