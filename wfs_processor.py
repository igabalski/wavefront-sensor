import warnings
warnings.filterwarnings("ignore")
import numpy as np
from scipy.ndimage.measurements import center_of_mass
from scipy.optimize import curve_fit
from scipy.ndimage import maximum_filter
from scipy.sparse.linalg import lsmr as iterative_solver
from scipy.sparse.linalg import LinearOperator
from itertools import combinations

import matplotlib.pyplot as plt

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
    local_maxima[:,:int((aoi_size+1)/2)] = 0
    local_maxima[:,-int((aoi_size+1)/2):] = 0
    local_maxima[:int((aoi_size+1)/2),:] = 0
    local_maxima[-int((aoi_size+1)/2):,:] = 0

    
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
    # aoi_size = size of aois in pixels
    # pixel_size = size of pixels in meters
    # magnification = magnification as a raw number
    #Returns:
    # r0 = Fried parameter measured from fit
    
    ssf_array = np.array(ssf_array)
    d_sub = aoi_size*pixel_size*magnification
    def functional_form_full(r, r0):
        val = np.multiply(6.88/(d_sub**2)*(d_sub/r0)**(5/3)*1.32, np.add(np.divide(np.square(r),np.add(5.4,np.square(r))),np.multiply(0.557, np.power(r, 0.213))))
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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def infer_aoi_size(aoi_locations):
    #Inputs:
    # aoi_locations = unsorted aoi locations; Either a dict of the form {'x,y': [x_centroid, y_centroid]'}, or a numpy array of the form [x,y]
    #Returns:
    # aoi_size_x, aoi_size_y = floats indicating inferred aoi separations in x and y directions

    if(isinstance(aoi_locations, dict)):
        aoi_locations = np.array([[int(x) for x in aoi.split(',')] for aoi in aoi_locations])


    x_locations, y_locations = aoi_locations[:,0], aoi_locations[:,1]
    xmin, xmax, ymin, ymax = np.amin(x_locations), np.amax(x_locations), np.amin(y_locations), np.amax(y_locations)

    x_hist, x_bin_edges = np.histogram(x_locations, bins=int(xmax-xmin))
    y_hist, y_bin_edges = np.histogram(y_locations, bins=int(ymax-ymin))
    
    x_hist = x_hist-np.mean(x_hist)
    y_hist = y_hist-np.mean(y_hist)

    x_fft = np.fft.fft(x_hist)
    y_fft = np.fft.fft(y_hist)

    x_freq = np.fft.fftfreq(len(x_hist), d=1)
    y_freq = np.fft.fftfreq(len(y_hist), d=1)


    aoi_size_x = 1/np.abs(x_freq[np.argmax(x_fft)])
    aoi_size_y = 1/np.abs(y_freq[np.argmax(y_fft)])
    
    return aoi_size_x, aoi_size_y



def sort_aois(aoi_locations, aoi_size, buffer_size=3):
    #Inputs:
    # aoi_locations = unsorted aoi locations; Either a dict of the form {'x,y': [x_centroid, y_centroid]'}, or a numpy array of the form [x,y]
    # aoi_size = the size of a wavefront sensor subaperture in pixels ; If single value, represents square aoi. If tuple, represents (aoi_size_x, aoi_size_y)
    # buffer_size = pixel buffer in each direction to allow for slightly irregular spacing of aois
    #Returns:
    # sorted_aoi_locations = numpy array of sorted aoi_locations list; array iterates top left to bottom right, row by row
    
    if(isinstance(aoi_locations, dict)):
        aoi_locations = np.array([[int(x) for x in aoi.split(',')] for aoi in aoi_locations])



    if(isinstance(aoi_size, tuple)):
        aoi_size_x, aoi_size_y = aoi_size
    else:
        aoi_size_x = aoi_size
        aoi_size_y = aoi_size

    ymin, ymax = np.amin(aoi_locations[:,1]), np.amax(aoi_locations[:,1])
    
    sorted_list = []
    for y in np.arange(ymin, ymax+1, aoi_size_y):
        row = np.array(aoi_locations[np.logical_and(aoi_locations[:,1]+buffer_size>=y, aoi_locations[:,1]-buffer_size<=y)])
        sorted_row = row[np.argsort(row[:,0], axis=0)]
        sorted_list.append(sorted_row)
    
    sorted_aoi_locations = np.array([xy_tuple for row in sorted_list for xy_tuple in row])

    return sorted_aoi_locations


def get_current_aoi(aoi, aoi_locations):
    #Inputs:
    # aoi = the current aoi location of the form [x,y]
    # aoi_locations = numpy array of all aoi locations
    #Returns:
    # mask = numpy array mask indicating which index in aoi_locations represents current aoi
    
    mask = [np.array_equal(x, aoi) for x in aoi_locations]
    
    return mask


def get_adjacent_aois(aoi, aoi_locations, aoi_size, buffer_size=3):
    #Inputs:
    # aoi = numpy 1-D array of the form [x,y] representing aoi pixel top left corner location
    # aoi_locations = numpy array of all aoi locations of the form [x,y]
    # aoi_size = the size of a wavefront sensor subaperture in pixels (assumes square AOIs)
    # buffer_size = pixel buffer in each direction to allow for slightly irregular spacing of aois
    #Returns:
    # mask = numpy array mask indicating which aois in aoi_locations are adjacent to aoi
    
    x_mask_lower = np.logical_and(aoi_locations[:,0]+aoi_size+buffer_size>=aoi[0], aoi_locations[:,0]+aoi_size-buffer_size<=aoi[0])
    x_mask_upper = np.logical_and(aoi_locations[:,0]-aoi_size+buffer_size>=aoi[0], aoi_locations[:,0]-aoi_size-buffer_size<=aoi[0])
    x_mask_equal = np.logical_and(aoi_locations[:,0]+buffer_size>=aoi[0], aoi_locations[:,0]-buffer_size<=aoi[0])
    y_mask_lower = np.logical_and(aoi_locations[:,1]+aoi_size+buffer_size>=aoi[1], aoi_locations[:,1]+aoi_size-buffer_size<=aoi[1])
    y_mask_upper = np.logical_and(aoi_locations[:,1]-aoi_size+buffer_size>=aoi[1], aoi_locations[:,1]-aoi_size-buffer_size<=aoi[1])
    y_mask_equal = np.logical_and(aoi_locations[:,1]+buffer_size>=aoi[1], aoi_locations[:,1]-buffer_size<=aoi[1])
    
    x_mask = np.logical_and(np.logical_or(x_mask_lower, x_mask_upper), y_mask_equal)
    y_mask = np.logical_and(np.logical_or(y_mask_lower, y_mask_upper), x_mask_equal)
    mask = np.logical_or(x_mask, y_mask)
    
    return mask


def get_aoi_signature(aoi, adjacent_aoi, buffer_size=3):
    #Inputs:
    # aoi = current aoi location of the form [x,y]
    # adjacent_aoi = adjacent aoi location to be compared to current aoi
    # buffer_size = pixel buffer in each direction to allow for slightly irregular spacing of aois
    #Returns:
    # signature = tuple of the form (x_signature, y_signature)
    # NOTE: each signature is +1 if adjacent_aoi is positively located w.r.t. aoi,
    #       -1 if negatively located w.r.t. aoi, else 0
    
    x_signature, y_signature = (0,0)
    
    if(adjacent_aoi[0]-buffer_size>aoi[0]+buffer_size):
        x_signature = 1
    elif(adjacent_aoi[0]+buffer_size<aoi[0]-buffer_size):
        x_signature = -1
        
    if(adjacent_aoi[1]-buffer_size>aoi[1]+buffer_size):
        y_signature = 1
    elif(adjacent_aoi[1]+buffer_size<aoi[1]-buffer_size):
        y_signature = -1
    
    signature = (x_signature, y_signature)
    
    return signature


def get_aoi_index(aoi, xmin, ymin, aoi_size, buffer_size=3):
    #Inputs:
    # aoi = current aoi location of the form [x,y]
    # xmin, ymin = minimum x and y aoi locations (does not have to be from same aoi)
    # aoi_size = the size of a wavefront sensor subaperture in pixels (assumes square AOIs)
    # buffer_size = pixel buffer in each direction to allow for slightly irregular spacing of aois
    #Returns:
    # index = array of the form [x_index, y_index] representing aoi location index relative to top left
    
    index = (int((aoi[0]-xmin+2*buffer_size)/aoi_size), int((aoi[1]-ymin+2*buffer_size)/aoi_size))
    
    return index


def to_string(aoi):
    #Inputs:
    # aoi = numpy int32 array of the form [x,y]
    #Returns:
    # string_representation = hashable string representation of aoi; use as key for gradients dict
    
    string_representation = str(aoi[0])+','+str(aoi[1])
    return string_representation


def build_reconstruction_matrix(sorted_aoi_locations, aoi_size):
    #Inputs:
    # sorted_aoi_locations = numpy array of sorted aoi_locations list; array iterates top left to bottom right, row by row
    # aoi_size = the size of a wavefront sensor subaperture in pixels (assumes square AOIs)
    #Returns:
    # A = Southwell reconstruction matrix, shape is (num_equations, num_aois)
    
    num_aois = len(sorted_aoi_locations)
    num_equations = num_aois+1
    
    A = np.zeros((num_equations, num_aois), dtype=np.float64)
    A[-1] = np.ones((num_aois))
    for i, aoi in enumerate(sorted_aoi_locations):
        adjacent_aois = sorted_aoi_locations[get_adjacent_aois(aoi, sorted_aoi_locations, aoi_size)]
        num_adjacent_aois = len(adjacent_aois)
        A[i, get_current_aoi(aoi, sorted_aoi_locations)]=-1
        if(num_adjacent_aois!=0):
            A[i, get_adjacent_aois(aoi, sorted_aoi_locations, aoi_size)]=1/num_adjacent_aois

    return A


def build_slope_vector(gradients, sorted_aoi_locations, aoi_size, pixel_size, magnification):
    #Inputs:
    # gradients = dictionary of form {'x_location,y_location': [x_gradient, y_gradient]} referenced from top left 
    # sorted_aoi_locations = numpy array of sorted aoi_locations list; array iterates top left to bottom right, row by row
    # aoi_size = the size of a wavefront sensor subaperture in pixels (assumes square AOIs)
    # pixel_size = size of pixels in meters
    # magnification = magnification as a raw number
    #Returns:
    # S = Southwell slope vector, shape is (num_equations)
    
    num_aois = len(sorted_aoi_locations)
    num_equations = num_aois+1
    aoi_separation = aoi_size*pixel_size*magnification
    
    S = np.zeros((num_equations))
    for i, aoi in enumerate(sorted_aoi_locations):
        adjacent_aois = sorted_aoi_locations[get_adjacent_aois(aoi, sorted_aoi_locations, aoi_size)]
        num_adjacent_aois = len(adjacent_aois)
        
        if(num_adjacent_aois!=0):
            current_gradient = gradients[to_string(aoi)]
            for adjacent_aoi in adjacent_aois:
                adjacent_gradient = gradients[to_string(adjacent_aoi)]
                signature = get_aoi_signature(aoi, adjacent_aoi)
                if(signature[1]==0):
                    avg_slope_val = signature[0]*aoi_separation*(adjacent_gradient[0]+current_gradient[0])/2
                elif(signature[0]==0):
                    avg_slope_val = signature[1]*aoi_separation*(adjacent_gradient[1]+current_gradient[1])/2
                S[i] += avg_slope_val
            S[i] /= num_adjacent_aois
    
    return S


def reconstruct_wavefront(A, S, sorted_aoi_locations, aoi_size, buffer_size=3):
    #Inputs:
    # A = Southwell reconstruction matrix, shape is (num_equations, num_aois)
    # S = Southwell slope vector, shape is (num_equations)
    # sorted_aoi_locations = numpy array of sorted aoi_locations list; array iterates top left to bottom right, row by row
    # aoi_size = the size of a wavefront sensor subaperture in pixels (assumes square AOIs)
    # buffer_size = pixel buffer in each direction to allow for slightly irregular spacing of aois
    #Returns:
    # wavefront = a numpy float32 array of the reconstructed wavefront (value is 0 if no aoi at location)
    # NOTE: Solves the system A*phi=S through iterative least-squares. Phi vector is then reshaped to 2-D wavefront image.
    
    phi, info, num_iterations, normr, normar, norma, conda, normx = iterative_solver(A, S)
    
    xmin, xmax = np.amin(sorted_aoi_locations[:,0]), np.amax(sorted_aoi_locations[:,0])
    ymin, ymax = np.amin(sorted_aoi_locations[:,1]), np.amax(sorted_aoi_locations[:,1])
    
    num_x, num_y = int((xmax-xmin+2*buffer_size)/aoi_size)+1, int((ymax-ymin+2*buffer_size)/aoi_size)+1
    
    wavefront = np.zeros((num_y, num_x))
    for i, aoi in enumerate(sorted_aoi_locations):
        index = get_aoi_index(aoi, xmin, ymin, aoi_size)
        wavefront[index[1], index[0]]=phi[i]
    
    return wavefront


def process_file(filepath, aoi_size, focal_length, pixel_size, wavelength, magnification, calculate_turbulence=False, reconstruct=False):
    #Inputs:
    # filepath = path to the file to be processed
    # aoi_size = the size of a wavefront sensor subaperture in pixels (assumes square AOIs)
    # focal_length = float representing focal length of lenslets in meters
    # pixel_size = size of pixels in meters
    # wavelength = wavelength in meters
    # magnification = magnification of telescope system as a raw number
    # calculate_turbulence = whether or not to calculate turbulence parameters with slope structure function methods
    # reconstruct = whether or not to reconstruct wavefront using Southwell reconstructor
    #Yields:
    # (Always) status = one of 'aoi locations', 'framenum', 'turbulence', 'wavefront', 'both', indicates what is being returned at each yield
    # (First) aoi_locations = a list of aoi location strings of the form 'x_topleft_corner,y_topleft_corner'
    # (Intermediate) framenum = index of frame that was just processed
    # (Final) turbulence_outputs = list of turbulence parameters and slope structure functions (possibly empty if turbulence not calculated)
    # (Final) wavefronts_list = list of reconstructed wavefront images (possibly empty if wavefront not reconstructed)

    '''
    NOTE:
    This function is a generator which iterates over frames and yields in between frames.
    The first yield value is the list of aoi location coordinates, referenced to the top left corner of the aoi.
    The intermediate yield values are the current frame index.
    The final yield value is a tuple of the final outputs.

    Proper use is:
    file_processor = process_file(filepath, aoi_size, focal_length, pixel_size, wavelength, magnification, calculate_turbulence=ct, reconstruct=r)
    for frame in file_processor:
        status, return_values = frame
        if(status=='aoi locations'):
            aoi_locations = np.array([[int(val) for val in line.split(',')] for line in return_values])
            ...update gui with aoi locations...
        elif(status=='framenum'):
            ...update gui with frame number...
        else:
            ...update gui with processed information...
    turbulence_outputs, wavefronts = return_values
    '''

    images_array, time_list, version, bitdepth = read_file(filepath)
    aoi_locations, summed_array = get_aoi_locations(images_array, aoi_size)
    status = 'aoi locations'
    yield status, aoi_locations
    references = calculate_references(summed_array, aoi_locations, aoi_size)
    ssf_list = []
    ssfx_list = []
    ssfy_list = []
    turbulence_outputs = []
    wavefronts_list = []
    aoi_size_tuple = infer_aoi_size(references)
    sorted_aoi_locations = sort_aois(references, aoi_size_tuple)
    A = build_reconstruction_matrix(sorted_aoi_locations, aoi_size)

    try:
        framenum = 0
        for frame in images_array:
            centroids = find_centroids(frame, aoi_locations, aoi_size)
            differences, gradients = find_slopes(centroids, references, focal_length, pixel_size, wavelength, magnification)
            if(calculate_turbulence):
                unbinned_ssf, unbinned_ssfx, unbinned_ssfy = get_unbinned_ssf(gradients, aoi_size)
                ssf = get_binned_ssf(unbinned_ssf)
                ssfx = get_binned_ssf(unbinned_ssfx)
                ssfy = get_binned_ssf(unbinned_ssfy)
                ssf_list.append(ssf)
                ssfx_list.append(ssfx)
                ssfy_list.append(ssfy)
            if(reconstruct):
                S = build_slope_vector(gradients, sorted_aoi_locations, aoi_size, pixel_size, magnification)
                wavefront = reconstruct_wavefront(A, S, sorted_aoi_locations, aoi_size)
                wavefronts_list.append(wavefront)
            status = 'framenum'
            yield status, framenum
            framenum += 1

        return_values = []
        
        if(calculate_turbulence):
            ssf_list, ssfx_list, ssfy_list = np.array(ssf_list), np.array(ssfx_list), np.array(ssfy_list)
            ssf = np.mean(ssf_list, axis=0)
            ssfx = np.mean(ssfx_list, axis=0)
            ssfy = np.mean(ssfy_list, axis=0)
            
            r0_full = fit_ssf(ssf, aoi_size, pixel_size, magnification, fit_type='full')
            r0_individual = np.mean((fit_ssf(ssfx, aoi_size, pixel_size, magnification, fit_type='individual'), fit_ssf(ssfy, aoi_size, pixel_size, magnification, fit_type='individual')))
            
            turbulence_outputs = (r0_full, r0_individual, ssf, ssfx, ssfy)
            status = 'turbulence'
        
        if(reconstruct):
            status = 'wavefront'

        if(calculate_turbulence and reconstruct):
            status = 'both'

        return_values.append(turbulence_outputs)
        return_values.append(wavefronts_list)

        yield status, return_values
    
    except StopIteration:
        pass



# CURRENT STATUS: implemented aoi size inferment algorithm
def run():
    filepath = '/home/ian/Desktop/workspace/example_acs_processing/data2018_08_28_16_11_29/data2018_08_28_16_11_29.dat'
    aoi_size = 20
    focal_length = 6.7e-3
    pixel_size = 7.4e-6
    wavelength = 640e-9
    magnification = 40
    calculate_turbulence=False
    reconstruct=True

    file_processor = process_file(filepath, aoi_size, focal_length, pixel_size, wavelength, magnification, calculate_turbulence=calculate_turbulence, reconstruct=reconstruct)
    
    status, return_values = next(file_processor)
    if(status=='aoi locations'):
        aoi_locations = np.array([[int(val) for val in line.split(',')] for line in return_values])
        aoi_size = infer_aoi_size(aoi_locations)
        sorted_aoi_locations = sort_aois(aoi_locations, aoi_size)




if __name__=='__main__':
    run()