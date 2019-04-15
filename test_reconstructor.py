import numpy as np
from numpy.random import randint
import wfs_processor as wfs
import matplotlib.pyplot as plt


aoi_size = 20
pixel_size = 1
magnification = 1

aoi_locations = np.array([np.array([(aoi_size*x), (aoi_size*y)]) for x in range(20) for y in range(20)])
randarray = randint(-1, 2, size=aoi_locations.shape)
aoi_locations = np.add(aoi_locations, randarray)

def gaussian_gradient(location, mu, sigma):
    xs, ys = location[0]-mu[0], location[1]-mu[1]
    return (-2/sigma**2*xs*np.exp(-(xs**2+ys**2)/sigma**2), -2/sigma**2*ys*np.exp(-(xs**2+ys**2)/sigma**2))

xmin, xmax = np.amin(aoi_locations[:,0]), np.amax(aoi_locations[:,0])
ymin, ymax = np.amin(aoi_locations[:,1]), np.amax(aoi_locations[:,1])
gradients={}
for aoi in aoi_locations:
    #print(aoi)
    aoi_string = str(aoi[0])+','+str(aoi[1])
    gradients[aoi_string] = gaussian_gradient(wfs.get_aoi_index(aoi, xmin, ymin, aoi_size), (10,10), 5)

aoi_locations = wfs.sort_aois(gradients, aoi_size)


sorted_aoi_locations = wfs.sort_aois(gradients, aoi_size, buffer_size=2)
A = wfs.build_reconstruction_matrix(sorted_aoi_locations, aoi_size)
S = wfs.build_slope_vector(gradients, sorted_aoi_locations, aoi_size, pixel_size, magnification)
wavefront = wfs.reconstruct_wavefront(A, S, sorted_aoi_locations, aoi_size)

plt.figure()
plt.imshow(wavefront)
plt.show()

    


