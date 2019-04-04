import warnings
warnings.filterwarnings("ignore")
import numpy as np
from scipy.ndimage.measurements import center_of_mass
from scipy.optimize import curve_fit
from scipy.ndimage import maximum_filter
from itertools import combinations


def read_file(filepath):
    #Inputs:
    # filepath = path to wavefront sensor image file to be read
    #Returns:
    # images_array = numpy array of image frames, 16-bit unsigned integer pixel values
    # time_list = time stamps for each frame
    # version = file format version
    # bitdepth = camera bit depth

    try:
        
        with open(filepath,"rb") as f:

            version = np.fromfile(f,count=1,dtype=np.double)
            bitdepth = np.fromfile(f,count=1,dtype=np.int32)
            
            time_list = []
            images_array = []

            while(True):
                frameID = np.fromfile(f,count=1,dtype=np.uint64)
                if(len(frameID)==0):
                    break
                time = np.fromfile(f,count=1,dtype=np.double)
                time_list.append(time)
                width = int(np.fromfile(f,count=1,dtype=np.uint32))
                height = int(np.fromfile(f,count=1,dtype=np.uint32))
                buffer_size = np.fromfile(f,count=1,dtype=np.uint32)
                #buffer size is size of a single frame in bytes. Frame reads out as a stream of pixel values.
                
                image_data = np.fromfile(f,count=int(width*height),dtype=np.uint16)
                if(len(image_data) == (width*height)):
                    image = image_data.reshape([height,width])
                    images_array.append(image)
                else:
                    print("Frame %s Skipped!!!" % (frameID))
        

        return np.array(images_array), np.array(time_list), version, bitdepth
            
    except Exception as e:
        print(e)
        return None


def get_aoi_locations(images_array, aoi_size, threshold=0.2):
    #Inputs:
    # images_array = the numpy array of images for which to calculate reference spots
    # aoi_size = the size of a wavefront sensor subaperture in pixels (assumes square AOIs)
    # threshold = all pixels in summed array less than maxval*threshold are set to 0
    #Returns:
    # aoi_locations = list of strings of the form 'x,y' where x, y are integer locations of top left corner of aois
    # summed_array = summed array, summed along each pixel axis individually
    
    summed_array = np.sum(images_array, axis=0)
    maxval = summed_array.max()
    pixel_threshold = maxval*threshold
    summed_array[summed_array<pixel_threshold] = 0
    # find highest regional maxima in aoi_sized window
    # this relies on there being exactly one true maximum in any aoi-sized window
    # if two pixels in summed_array within a window are exactly the same value, they will both be treated as aoi locations (this almost never happens)

    maxfiltered_image=maximum_filter(summed_array, size=aoi_size)
    maxfiltered_image_reversed = maximum_filter(summed_array[::-1,::-1], size=aoi_size)
    local_maxima = np.array(((summed_array==maxfiltered_image)*(summed_array==maxfiltered_image_reversed[::-1,::-1])*(summed_array>0)), dtype=np.int32)
    local_maxima[:,:int((aoi_size+1)/2)]=0
    local_maxima[:,-int((aoi_size+1)/2):]=0
    local_maxima[:int((aoi_size+1)/2),:]=0
    local_maxima[-int((aoi_size+1)/2):,:]=0

    
    (xcoords, ycoords) = np.where(local_maxima==1)
    local_maxima_coords = np.transpose(np.stack((ycoords, xcoords)))

    #subtract half an aoi size from local maxima locations to get top left corners of aois
    aoi_offset = int((aoi_size-1)/2)
    aoi_locations_array = local_maxima_coords - np.full(shape=local_maxima_coords.shape, fill_value=aoi_offset)

    aoi_locations = [str(x[0])+','+str(x[1]) for x in aoi_locations_array]

    return aoi_locations, summed_array


def find_centroids(image, aoi_locations, aoi_size):
    #Inputs:
    # image = numpy ndarray of floating point values
    # aoi_locations = list of strings of the form 'x,y' where x, y are integer locations of top left corner of aois
    # aoi_size_x, aoi_size_y = aoi size in pixels
    #Returns:
    # centroids = dictionary of the form {'x_location,y_location':[x_centroid, y_centroid]} referenced from top left of aois
    
    centroids={}
    for aoi in aoi_locations:
        xmin=int(aoi.split(',')[0])
        ymin=int(aoi.split(',')[1])
        xmax, ymax = xmin+aoi_size, ymin+aoi_size
        centroid=center_of_mass(image[ymin:ymax, xmin:xmax])
        centroids[aoi]=centroid

    return centroids


def calculate_references(summed_array, aoi_locations, aoi_size):
    #Inputs:
    # summed_array = the summed image of wavefront sensor spots from which to calculate reference centroids
    # aoi_locations = list of strings of the form 'x,y' where x, y are integer locations of top left corner of aois
    # aoi_size = the size of a wavefront sensor subaperture in pixels (assumes square AOIs)
    # threshold = all pixels less than the max pixel value times threshold set to zero to eliminate noise
    #Returns:
    # references = a dictionary of the form {'x_location,y_location':[x_centroid, y_centroid]} referenced from top left of aois, representing reference centroids

    references = find_centroids(summed_array, aoi_locations, aoi_size)


    return references


def find_slopes(centroids, references, focal_length, pixel_size, wavelength, magnification):
    #Inputs:
    # centroids = dictionary of form {'x_location,y_location': [x_centroid, y_centroid]} referenced from top left
    # references = dictionary of form {'x_location,y_location': [x_reference, y_reference]} referenced from top left
    # focal_length = float representing focal length of lenslets in meters
    # wavelength = wavelength in meters
    # magnification = magnification as a raw number
    # pixel_size = float representing pixel pitch in meters
    #Returns:
    # differences = dictionary of the form {'x_location,y_location':[relative_x_centroid, relative_y_centroid]} in pixels
    # gradients = dictionary of the form {'x_location,y_location':[phase_x_gradient, phase_y_gradient]} in meters^-1
    
    differences = {}
    gradients = {}
    for location in references.keys():
        differences[location] = np.subtract(centroids[location], references[location])
        factor = 2*np.pi/wavelength/focal_length/magnification*pixel_size
        gradients[location] = np.multiply(factor, differences[location])
    
    return differences, gradients


def get_unbinned_ssf(gradients, aoi_size):
    #Inputs:
    # gradients = dictionary of form {'x_location,y_location': [x_gradient, y_gradient]} referenced from top left
    # aoi_size = the size of a wavefront sensor subaperture in pixels (assumes square AOIs)
    #Returns:
    # ssf = unbinned full slope structure function array of the form [r_separation (aois), (grad_phi1-grad_phi2)**2]
    # ssfx = unbinned x slope structure function array of the same form
    # ssfy = unbinned y slope structure function array of the same form
    
    ssf = []
    ssfx = []
    ssfy = []
    i=0
    
    gradients_list = list(gradients.items())
    gradients_array = np.array([[[float(x) for x in gradient[0].split(',')], gradient[1]] for gradient in gradients_list])
    aoi_pairs = np.array(list(combinations(gradients_array, 2)))
    
    pair_difference = np.subtract(aoi_pairs[:,0], aoi_pairs[:,1])
    separations = np.linalg.norm(np.multiply(1/aoi_size, np.subtract(aoi_pairs[:,0,0], aoi_pairs[:,1,0])), axis=1)
    xpair_difference_squared = np.square(pair_difference[:,1,0])
    ypair_difference_squared = np.square(pair_difference[:,1,1])
    pair_difference_squared = xpair_difference_squared+ypair_difference_squared
    
    ssf = np.transpose(np.stack((separations, pair_difference_squared)))
    ssfx = np.transpose(np.stack((separations, xpair_difference_squared)))
    ssfy = np.transpose(np.stack((separations, ypair_difference_squared)))

    
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


def fit_ssf(ssf_array, aoi_size, pixel_size, magnification, fit_type='full'):
    #Inputs:
    # ssf_array = slope structure function array of the form [r_separation (aois), (grad_phi1-grad_phi2)**2]
    # aoi_size_x, aoi_size_y = x and y sizes of aois in pixels
    # pixel_size = size of pixels in meters
    # magnification = magnification as a raw number
    #Returns:
    # r0 = Fried parameter measured from fit
    
    ssf_array = np.array(ssf_array)
    d_sub = aoi_size*pixel_size*magnification
    def functional_form_full(r, r0):
        val = np.multiply(6.88/(d_sub**2)*(d_sub/r0)**(5/3)*1.32, np.add(np.divide(np.square(r),np.add(5.4,np.square(r))),np.multiply(0.557, np.power(r, 0.213))))
        #print(val)
        return val

    def functional_form_individual(r, r0):
        val = np.multiply(6.88/(d_sub**2)*(d_sub/r0)**(5/3)*0.77, np.add(np.divide(np.square(r),np.add(1,np.square(r))),np.multiply(0.438, np.power(r, 0.2555))))
        return val

    if(fit_type=='individual'):
        functional_form = functional_form_individual
    else:
        functional_form = functional_form_full

        
    params, covariances = curve_fit(functional_form, ssf_array[:,0], ssf_array[:,1], p0=0.01, bounds=(0.0001, 1))
    return params[0]


def process_all_frames(images_array, references, aoi_size, focal_length, pixel_size, wavelength, magnification):
    #Inputs:
    # images_array = the numpy array of images for which to calculate reference spots
    # reference_spots = a dictionary of the form {'x_location,y_location':[x_centroid, y_centroid]} referenced from top left of aois, representing reference centroids
    # aoi_size = the size of a wavefront sensor subaperture in pixels (assumes square AOIs)
    # focal_length = float representing focal length of lenslets in meters
    # pixel_size = size of pixels in meters
    # wavelength = wavelength in meters
    # magnification = magnification of telescope system as a raw number
    #Returns:
    # r0_full = Fried parameter r0 calculated using full slope structure function method
    # r0_individual = Fried parameter calculated using individual slope structure functions separately and averaging results
    # ssf = full slope structure function numpy array of the form [r, ssf(r)] where r is measured in aois
    # ssfx = x slope structure function numpy array of the form [r, ssfx(r)] where r is measured in aois
    # ssfy = fy slope structure function numpy array of the form [r, ssfy(r)] where r is measured in aois

    ssf_list = []
    ssfx_list = []
    ssfy_list = []


    aoi_locations = list(references.keys())
    for frame in images_array:
        centroids = find_centroids(frame, aoi_locations, aoi_size)
        differences, gradients = find_slopes(centroids, references, focal_length, pixel_size, wavelength, magnification)
        unbinned_ssf, unbinned_ssfx, unbinned_ssfy = get_unbinned_ssf(gradients, aoi_size)
        ssf = get_binned_ssf(unbinned_ssf)
        ssfx = get_binned_ssf(unbinned_ssfx)
        ssfy = get_binned_ssf(unbinned_ssfy)
        ssf_list.append(ssf)
        ssfx_list.append(ssfx)
        ssfy_list.append(ssfy)

    ssf_list, ssfx_list, ssfy_list = np.array(ssf_list), np.array(ssfx_list), np.array(ssfy_list)
    ssf = np.mean(ssf_list, axis=0)
    ssfx = np.mean(ssfx_list, axis=0)
    ssfy = np.mean(ssfy_list, axis=0)
    

    r0_full = fit_ssf(ssf, aoi_size, pixel_size, magnification, fit_type='full')
    r0_individual = np.mean((fit_ssf(ssfx, aoi_size, pixel_size, magnification, fit_type='individual'), fit_ssf(ssfy, aoi_size, pixel_size, magnification, fit_type='individual')))

    return r0_full, r0_individual, ssf, ssfx, ssfy

def process_file(filepath, aoi_size, focal_length, pixel_size, wavelength, magnification):
    #Inputs:
    # filepath = path to wavefront sensor raw video file to be processed
    # aoi_size = the size of a wavefront sensor subaperture in pixels (assumes square AOIs)
    # focal_length = float representing focal length of lenslets in meters
    # pixel_size = size of pixels in meters
    # wavelength = wavelength in meters
    # magnification = magnification of telescope system as a raw number
    #Returns:
    # r0_full = Fried parameter r0 calculated using full slope structure function method
    # r0_individual = Fried parameter calculated using individual slope structure functions separately and averaging results
    # ssf = full slope structure function numpy array of the form [r, ssf(r)] where r is measured in aois
    # ssfx = x slope structure function numpy array of the form [r, ssfx(r)] where r is measured in aois
    # ssfy = y slope structure function numpy array of the form [r, ssfy(r)] where r is measured in aois

    images_array, time_list, version, bitdepth = read_file(filepath)
    aoi_locations, summed_array = get_aoi_locations(images_array, aoi_size)
    references = calculate_references(summed_array, aoi_locations, aoi_size)
    r0_full, r0_individual, ssf, ssfx, ssfy = process_all_frames(images_array, references, aoi_size, focal_length, pixel_size, wavelength, magnification)

    return r0_full, r0_individual, ssf, ssfx, ssfy