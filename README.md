# Divisive-Normalization-Transform
implementation of Divisive Normalization Transform (DNT) with Python.
Divisive normalization has been recognized as a successful approach to model the perceptualsensitivity of biological vision. It also provides a useful image representation that significantly improves statistical independence for natural images. The divisive normalization procedure tends to produce coefficients distributed  in  a  Gaussian  manner  for  natural  images.  In  the presence  of  perturbation  however,  this  Gaussianity  at  the output  of  the  normalization  procedure  is  not  guaranteed. 

The algorithm is described in:
Li, Qiang and Wang, Zhou, "Reduced-reference image quality assessment using divisive normalization-based image representation", IEEE journal of selected topics in signal processing.

# Parameters:
After having the shape and the coeficients of subbands using the steerable pyramid wavelet-transform for each scale in each orientations.

pyro: coeficients of all the subband of the wavelet-transform.
pind: shape of each subband.
num_scales: number of the used scales.
num_or: number of orientation in each scale.

# Returns:
subband: num_scales*num_or subband.
shape: a new shape for each subband.  
