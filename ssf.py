import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.ndimage.measurements import center_of_mass
from scipy.optimize import curve_fit
from itertools import combinations

def read_info_file(filename):
	#Inputs:
	# filename = name of information file, expects a .txt file
	# file is of format:
	# aoi_size_x (pixels)
	# aoi_size_y (pixels)
	# focal_length (meters)
	# pixel_size (meters)
	# magnification (number)
	# x_location, y_location (top left corner), x_reference, y_reference (reference relative to top left, float number of pixels)
	#Returns:
	# aoi_size_x = x aoi size in pixels
	# aoi_size_y = y aoi size in pixels
	# focal_length = focal length of lenslet array in meters
	# pixel_size = physical size of pixels in meters
	# wavelength = wavelength of light in meters
	# magnification = magnification of telescope-collimating lens system as a number
	# aois = dictionary of the form {'x_location,y_location':[x_reference, y_reference]} reference relative to top left
	with open(filename) as f:
		info = [x.strip() for x in f.read().split('\n')]
		aoi_size_x = int(info[0])
		aoi_size_y = int(info[1])
		focal_length = float(info[2])
		pixel_size = float(info[3])
		wavelength = float(info[4])
		magnification = float(info[5])
		pixel_size = 0.9766e-6
		aoi_size_x = 204.8
		aoi_size_y = 204.8
		magnification = 35
		
		aois={}
		for aoi_string in info[6:]:
			aoi=[x.strip() for x in aoi_string.split(',')]
			if(len(aoi)>1):
				location=aoi[0]+','+aoi[1]
				reference=[float(val) for val in aoi[2:4]]
				aois[location]=reference
	return aoi_size_x, aoi_size_y, focal_length, pixel_size, wavelength, magnification, aois

def read_image_file(filename):
	#Inputs:
	# filename = name of image file with formatting, expects a .mat file
	#Returns:
	# image['irradianceArray'] = image array object containing float-type pixel values
	image=loadmat(filename)
	return np.array(image[list(image.keys())[-1]])

def find_centroids(image, aoi_locations, aoi_size_x, aoi_size_y):
	#Inputs:
	# image = numpy ndarray of doubles
	# aoi_locations = list of strings of the form 'x,y' where x, y are integer locations of top left corner of aois
	# aoi_size_x, aoi_size_y = aoi size in pixels
	#Returns:
	# centroids = dictionary of the form {'x_location,y_location':[x_centroid, y_centroid]} referenced from top left
	
	centroids={}
	for aoi in aoi_locations:
		xmin=int(aoi.split(',')[0])
		ymin=int(aoi.split(',')[1])
		xmax, ymax = xmin+int(aoi_size_x)+1, ymin+int(aoi_size_y)+1
		centroid=np.array(center_of_mass(image[ymin:ymax, xmin:xmax]))
		centroids[aoi]=centroid
	return centroids

def find_slopes(centroids, references, focal_length, pixel_size, wavelength, magnification):
	#Inputs:
	# centroids = dictionary of form {'x_location,y_location': [x_centroid, y_centroid]} referenced from top left
	# references = dictionary of form {'x_location,y_location': [x_reference, y_reference]} referenced from top left
	# focal_length = float representing focal length in meters
	# wavelength = wavelength in meters
	# magnification = magnification as a raw number
	# pixel_size = float representing pixel pitch in meters
	#Returns:
	# differences = dictionary of the form {'x_location,y_location':[relative_x_centroid, relative_y_centroid]} in pixels
	# gradients = dictionary of the form {'x_location,y_location':[phase_x_gradient, phase_y_gradient]} in meters^-1
	
	differences={}
	gradients={}
	for location, centroid in centroids.items():
		differences[location]=np.subtract(centroid, references[location])
		factor=2*np.pi/wavelength/focal_length/magnification*pixel_size
		gradients[location]=np.multiply(factor, differences[location])
	
	return differences, gradients

def get_unbinned_ssf(gradients, aoi_size_x, aoi_size_y):
	#Inputs:
	# gradients = dictionary of form {'x_location,y_location': [x_gradient, y_gradient]} referenced from top left
	# aoi_size_x, aoi_size_y = aoi sizes in pixels
	#Returns:
	# ssf = unbinned full slope structure function array of the form [r_separation (aois), (grad_phi1-grad_phi2)**2]
	# ssfx = unbinned x slope structure function array of the same form
	# ssfy = unbinned y slope structure function array of the same form
	
	aoi_size=np.mean([aoi_size_x, aoi_size_y])
	ssf=[]
	ssfx=[]
	ssfy=[]
	for location_pair in combinations(gradients, 2):
		r1=[int(x) for x in location_pair[0].split(',')]
		r2=[int(x) for x in location_pair[1].split(',')]
		pair_difference=np.subtract(gradients[location_pair[0]], gradients[location_pair[1]])
		xpair_difference_squared=np.square(pair_difference[0])
		ypair_difference_squared=np.square(pair_difference[1])
		pair_difference_squared=xpair_difference_squared+ypair_difference_squared
		separation=np.linalg.norm(np.multiply(1/aoi_size, np.subtract(r1, r2)))
		ssf.append((separation, pair_difference_squared))
		ssfx.append((separation, xpair_difference_squared))
		ssfy.append((separation, ypair_difference_squared))
	
	return np.array(ssf), np.array(ssfx), np.array(ssfy)

def get_binned_ssf(unbinned_ssf, bin_size=1):
	#Inputs:
	# unbinned_ssf = numpy array, unbinned full slope structure function array of the form [r_separation (aois), (grad_phi1-grad_phi2)**2]
	# bin_size = size of bins in aois
	#Returns:
	# binned_ssf = binned slope structure function of the form [r_separation (aois), mean(grad_phi1-grad_phi2)**2(r_separation)]
	
	rmin = int(np.amin(unbinned_ssf[:,0]))
	rmax = int(np.amax(unbinned_ssf[:,0]))
	binned_ssf = [ [r+bin_size/2, np.mean(unbinned_ssf[np.logical_and(unbinned_ssf[:,0]>=r, unbinned_ssf[:,0]<r+bin_size)][:,1])] for r in range(rmin, rmax, bin_size)]
	return np.array(binned_ssf)

def fit_ssf(ssf_array, aoi_size_x, aoi_size_y, pixel_size, magnification):
	#Inputs:
	# ssf_array = slope structure function array of the form [r_separation (aois), (grad_phi1-grad_phi2)**2]
	# aoi_size_x, aoi_size_y = x and y sizes of aois in pixels
	# pixel_size = size of pixels in meters
	# magnification = magnification as a raw number
	#Returns:
	# r0 = Fried parameter measured from fit
	
	ssf_array=np.array(ssf_array)
	d_sub = np.mean([aoi_size_x, aoi_size_y])*pixel_size*magnification
	def functional_form(r, r0):
		val = np.multiply(6.88/(d_sub**2)*(d_sub/r0)**(5/3)*1.32, np.add(np.divide(np.square(r),np.add(5.4,np.square(r))),np.multiply(0.557, np.power(r, 0.213))))
		#print(val)
		return val
		
	params, covariances = curve_fit(functional_form, ssf_array[:,0], ssf_array[:,1], p0=0.01, bounds=(0.0001, 1))
	return params[0]
	
def fit_ssf_individual(ssf_array, aoi_size_x, aoi_size_y, pixel_size, magnification):
	#Inputs:
	# ssf_array = slope structure function array of the form [r_separation (aois), (grad_phi1-grad_phi2)**2]
	# aoi_size_x, aoi_size_y = x and y sizes of aois in pixels
	# pixel_size = size of pixels in meters
	# magnification = magnification as a raw number
	#Returns:
	# r0 = Fried parameter measured from fit
	
	ssf_array=np.array(ssf_array)
	d_sub = np.mean([aoi_size_x, aoi_size_y])*pixel_size*magnification
	def functional_form(r, r0):
		val = np.multiply(6.88/(d_sub**2)*(d_sub/r0)**(5/3)*0.77, np.add(np.divide(np.square(r),np.add(1,np.square(r))),np.multiply(0.438, np.power(r, 0.2555))))
		return val
		
	params, covariances = curve_fit(functional_form, ssf_array[:,0], ssf_array[:,1], p0=0.01, bounds=(0.0001, 1))
	return params[0]
	
def threshold_image(image, factor=0.2):
	#Inputs:
	# image = the image to be thresholded
	# factor = a number between 0 and 1 representing fraction of the image's max value to be thresholded out
	#Returns:
	# image = the thresholded image
	
	maxval=np.amax(image)
	image[image<factor*maxval]=0
	return image
	

def run():
	image = threshold_image(read_image_file('shwfsImageArray.mat'), factor=0)
	
	aoi_size_x, aoi_size_y, focal_length, pixel_size, wavelength, magnification, aois = read_info_file('shwfsImageData.txt')
	aoi_locations = list(aois.keys())
	aoi_references = aois
	
	centroids = find_centroids(image, aoi_locations, aoi_size_x, aoi_size_y)
	
	differences, gradients = find_slopes(centroids, aoi_references, focal_length, pixel_size, wavelength, magnification)
	
	unbinned_ssf, unbinned_ssfx, unbinned_ssfy = get_unbinned_ssf(gradients, aoi_size_x, aoi_size_y)
	binned_ssf = get_binned_ssf(unbinned_ssf)
	binned_ssfx = get_binned_ssf(unbinned_ssfx)
	binned_ssfy = get_binned_ssf(unbinned_ssfy)
		
	
	r0 = fit_ssf(binned_ssf, aoi_size_x, aoi_size_y, pixel_size, magnification)
	r0x = fit_ssf_individual(binned_ssfx, aoi_size_x, aoi_size_y, pixel_size, magnification)
	r0y = fit_ssf_individual(binned_ssfy, aoi_size_x, aoi_size_y, pixel_size, magnification)
	print(r0, np.mean([r0x, r0y]))
	x=np.arange(0,26,0.1)
	def functional_form(r, r0):
		d_sub=np.mean([aoi_size_x, aoi_size_y])*pixel_size*magnification
		val = np.multiply(6.88/(d_sub**2)*(d_sub/r0)**(5/3)*1.32, np.add(np.divide(np.square(r),np.add(5.4,np.square(r))),np.multiply(0.557, np.power(r, 0.213))))
		return val
	def functional_form2(r, r0):
		d_sub=np.mean([aoi_size_x, aoi_size_y])*pixel_size*magnification
		val = np.multiply(6.88/(d_sub**2)*(d_sub/r0)**(5/3)*0.77, np.add(np.divide(np.square(r),np.add(1,np.square(r))),np.multiply(0.438, np.power(r, 0.2555))))
		return val
	
	fitssf = functional_form(x, r0)
	fitssfx = functional_form2(x, r0x)
	fitssfy = functional_form2(x, r0y)
	theoryssf = functional_form(x, 0.0068202)
	theoryssf_indiv = functional_form2(x, 0.0068202)
	
	plt.figure()
	plt.plot(binned_ssf[:,0], binned_ssf[:,1], 'r--', label='Data')
	plt.plot(x, fitssf, 'r-', label='Fit')
	#plt.plot(x, theoryssf, 'b-', label='Theoretical')
	plt.legend()
	plt.savefig('xy_ssf.png')
	plt.show()
	
	plt.figure()
	plt.plot(binned_ssfx[:,0], binned_ssfx[:,1], 'r--', label='X Data')
	plt.plot(x, fitssfx, 'r-', label='X Fit')
	plt.plot(binned_ssfy[:,0], binned_ssfy[:,1], 'g--', label='Y Data')
	plt.plot(x, fitssfy, 'g-', label='Y Fit')
	#plt.plot(x, theoryssf_indiv, 'b-', label='Theoretical')
	plt.legend()
	plt.savefig('full_ssf.png')
	plt.show()
		
	


if __name__=='__main__':
	run()
# providing a comment here at the end for testing purposes
