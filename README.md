# Wavefront Sensor
This project provides Python code for post-processing a Shack-Hartmann wavefront sensor video file to calculate atmospheric turbulence parameters.
The code implements both the original slope structure function method, developed by Silbaugh et. al. (1996) and made computationally efficient by 
Mansell et. al. (2015), and the modified slope structure function method I developed recently. The two implemented 
methods share the same low-level Shack-Hartmann image processing methods and differ only by the particular slope structure function data they fit 
and the functional form used to perform the fit.
