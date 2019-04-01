import warnings
warnings.filterwarnings("ignore")
import numpy as np
from scipy.ndimage.measurements import center_of_mass
from scipy.optimize import curve_fit
from itertools import combinations

def calculate_reference_spots(self):
	self.status_label.config(text='Status: Calculating Reference Spots...')
	self.reference_spots=[]
	#sum and threshold all frames to average
	threshold=0.2
	self.summed_array_unthresholded=np.sum(self.imagesArray, axis=0)
	self.summed_array=np.array([i for i in self.summed_array_unthresholded])
	maxval=np.amax(self.summed_array)
	pixel_threshold=maxval*threshold
	self.summed_array[self.summed_array<pixel_threshold]=0
	
	
	self.im.set_data(self.summed_array)
	self.canvas.draw()
	test_window_size=int(self.controlframe.aoisize_entry.get())-5
	local_maxima=np.zeros((len(self.imagesArray[0]), len(self.imagesArray[0][0])))
	
	#find all local maxima in image
	for y in range(1, len(self.summed_array)-1):
		for x in range(1, len(self.summed_array[y])-1):
			val=self.summed_array[y][x]
			negx=self.summed_array[y][x-1]
			negy=self.summed_array[y-1][x]
			posx=self.summed_array[y][x+1]
			posy=self.summed_array[y+1][x]
			if(val>max(negx,negy,posx,posy)):
				local_maxima[y][x]=1
			if(not self.process_data):
				break
		if(not self.process_data):
			break
	
	#enforce constraint that there be only one local maximum within any 'test_window_size' (15x15) square
	#always take the highest local maximum when deciding
	for y in range(len(local_maxima)):
		for x in range(len(local_maxima[y])):
			if(local_maxima[y:y+test_window_size,x:x+test_window_size].sum()>1):
				maxval=self.summed_array[y][x]
				maxloc=[x,y]
				for i in range(y,y+test_window_size):
					for j in range(x,x+test_window_size):
						local_maxima[i][j]=0
						if(self.summed_array[i][j]>maxval):
							maxval=self.summed_array[i][j]
							maxloc=[j,i]
				local_maxima[maxloc[1]][maxloc[0]]=1
			if(not self.process_data):
				break
		if(not self.process_data):
			break
	for y in reversed(range(len(local_maxima)-test_window_size)):
		for x in reversed(range(len(local_maxima[y])-test_window_size)):
			if(local_maxima[y:y+test_window_size,x:x+test_window_size].sum()>1):
				maxval=self.summed_array[y][x]
				maxloc=[x,y]
				for i in range(y,y+test_window_size):
					for j in range(x,x+test_window_size):
						local_maxima[i][j]=0
						if(self.summed_array[i][j]>maxval):
							maxval=self.summed_array[i][j]
							maxloc=[j,i]
						if(not self.process_data):
							break
					if(not self.process_data):
						break
				local_maxima[maxloc[1]][maxloc[0]]=1
			if(not self.process_data):
				break
		if(not self.process_data):
			break
	
	ylim=len(local_maxima)
	xlim=len(local_maxima[0])
	
	local_maxima[:,0:int(test_window_size/2+5)]=0
	local_maxima[:,xlim-int(test_window_size/2+5):xlim]=0
	local_maxima[0:int(test_window_size/2+5),:]=0
	local_maxima[ylim-int(test_window_size/2+5):ylim,:]=0
	
	
	#add found reference spots to array, and red circles to image for testing
	window_size=int(self.controlframe.aoisize_entry.get())-3
	window_offset=int((window_size-1)/2) 
	for y in range(len(local_maxima)):
		for x in range(len(local_maxima[y])):
			if(local_maxima[y][x]==1):
				self.a.add_patch(patches.Circle((x,y),radius=2,edgecolor='r',facecolor='r'))
				self.a.add_patch(patches.Rectangle((x-window_offset-1,y-window_offset-1),window_size,window_size, edgecolor='c', facecolor='none'))
				
				rel_cent_x, rel_cent_y=center_of_mass(self.summed_array_unthresholded[y-window_offset:y-window_offset+window_size,x-window_offset:x-window_offset+window_size])
				cent_x=rel_cent_x+x-window_offset
				cent_y=rel_cent_y+y-window_offset
				
				self.reference_spots.append([x,y, cent_x, cent_y])
			if(not self.process_data):
				break
		if(not self.process_data):
			break
	self.canvas.draw()