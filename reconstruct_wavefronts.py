"""
Script to iterate over and process all wavefront sensor files in a directory.
Author: Ian Gabalski

This is a script which iterates over a directory tree and performs processing on
each wavefront sensor data file it finds in the tree. It also writes out a file
containing the names of the files which have already been processed, in case 
processing must be stopped and restarted.

"""

import wfs_processor as wfs
import os
import numpy as np
import struct


data_directory = '/home/ian/workspace/acs/data2018_08_28_16_09_18/'
wavefronts_directory = '/home/ian/workspace/acs/data2018_08_28_16_09_18/'
aoi_size = 20
focal_length = 6.7e-3
pixel_size = 7.4e-6
wavelength = 640e-9
magnification = 20


def write_wavefront_file(filepath, wavefronts):
	with open(filepath, 'wb') as filewriter:
		filewriter.write(struct.pack('d',3.0))
		filewriter.write(struct.pack('i', 16))
		for framenum, frame in enumerate(wavefronts):
			ymax, xmax = frame.shape
			filewriter.write(struct.pack('Q',framenum))				#frameID
			filewriter.write(struct.pack('d',0))					#timeIn
			filewriter.write(struct.pack('I',xmax))					#width
			filewriter.write(struct.pack('I',ymax))					#height
			filewriter.write(struct.pack('I',(xmax)*(ymax)*8))		#bufferSize
			for row in frame:
				for val in row:
					filewriter.write(struct.pack('d', val))


def process_files(filepath=data_directory, aoi_size=aoi_size, focal_length=focal_length, pixel_size=pixel_size, wavelength=wavelength, magnification=magnification):

	files = os.listdir(filepath)
	with open(data_directory+'processed_files.txt', 'w+')as f:
		f.write('Processed files:\n')

	for filenum, file in enumerate(files):
		if(file.startswith('data') and file.endswith('.dat')):
			print(file)
			print('File {} of {}: {}'.format(filenum, len(files), file))
			data_filepath = filepath+file
			file_processor = wfs.process_file(data_filepath, aoi_size, focal_length, pixel_size, wavelength, magnification, calculate_turbulence=False, reconstruct=True)
			for framenum, frame in enumerate(file_processor):
				status, return_values = frame
				if(status=='aoi locations'):
					aoi_locations = np.array([[int(val) for val in line.split(',')] for line in return_values])
					pass
				elif(status=='framenum'):
					pass
				else:
					pass
			turbulence_outputs, wavefronts = return_values
			with open(data_directory+'processed_files.txt', 'a+')as f:
				f.write(file+'\n')
			wavefront_filename = 'wavefront'+file[4:]
			wavefront_filepath = wavefronts_directory+wavefront_filename
			write_wavefront_file(wavefront_filepath, wavefronts)


if(__name__=='__main__'):
	process_files()
