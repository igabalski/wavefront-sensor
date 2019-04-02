import warnings
warnings.filterwarnings("ignore")
import numpy as np
from scipy.ndimage.measurements import center_of_mass
from scipy.optimize import curve_fit
from scipy.ndimage import maximum_filter
from itertools import combinations

def calculate_reference_spots(images_array, aoi_size, threshold=0.2):
	#Inputs:
	# images_array = the numpy array of images for which to calculate reference spots
	# aoi_size = the size of a wavefront sensor subaperture in pixels (assumes square AOIs)
	# threshold = all pixels less than the max pixel value times threshold set to zero to eliminate noise
	#Returns:
	# reference_spots = a list of reference spot integer and centroid locations of the form (x_int, y_int, x_cent, y_cent)
	reference_spots = []
	summed_array = np.sum(images_array, axis=0)
	maxval = np.amax(summed_array)
	pixel_threshold = maxval*threshold
	summed_array[summed_array<pixel_threshold]=0

	# find highest regional maxima in aoi_sized window
	# this relies on there being exactly one true maximum in any aoi-sized window
	# if two pixels in summed_array within a window are exactly the same value, they will both be treated as aoi locations (this almost never happens)

	maxfiltered_image=maximum_filter(summed_array, size=aoi_size, mode='nearest')
	maxfiltered_image_reversed = maximum_filter(summed_array[::-1,::-1], size=window_size, mode='nearest')
	local_maxima = np.array(((summed_array==maxfiltered_image)*(summed_array[::-1,::-1]==maxfiltered_image_reversed)[::-1,::-1]*(summed_array>0)), dtype=np.int32)
	
	(xcoords, ycoords) = np.where(local_maxima==1)
	local_maxima_coords = np.stack((xcoords, ycoords))

	for location in local_maxima_coords:
		xmin = np.amax((0,np.int(location[0]-(aoi_size/2))))
		ymin = np.amax((0,np.int(location[1]-(aoi_size/2))))
		xmax = np.amin((xmin+aoi_size,summed_array.shape[1]))
		ymax = np.amin((ymin+aoi_size,summed_array.shape[0]))
		x_centroid, y_centroid = center_of_mass(summed_array[ymin:ymax, xmin:xmax])
		reference_spots.append((location[0], location[1], x_centroid, y_centroid))

	return reference_spots