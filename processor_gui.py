import tkinter as tk
from tkinter.filedialog import askopenfilename, asksaveasfilename, askdirectory
import os
import time
import numpy as np
import threading

from tkinter import ttk
import matplotlib
import matplotlib.patches as patches
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib.backend_bases import key_press_handler
from mpl_toolkits.axes_grid1 import make_axes_locatable

import wfs_processor as wfs



class ACSDataApp(tk.Tk):

	def safe_quit(self):
		self.frames[ACSImageFrame].cancel_data_processing()
		time.sleep(0.2)
		self.quit()
		self.destroy()

	
	def __init__(self, *args, **kwargs):
		
		tk.Tk.__init__(self, *args, **kwargs)
		
		container=tk.Frame(self)
		container.pack(side='top', fill='both', expand=True)
		container.grid_rowconfigure(0, weight=1)
		container.grid_columnconfigure(0, weight=1)
		container.grid_columnconfigure(1, weight=15)
		
		self.wm_title('ACS Data Analysis')
		
		self.frames={}
		image_frame=ACSImageFrame(container, self)
		plots_frame=ACSPlotsFrame(container, self)
		controls_frame=ACSControlsFrame(container, self, image_frame, plots_frame)
		wavefront_frame=ACSWavefrontFrame(container, self)
		
		
		self.frames[ACSControlsFrame]=controls_frame
		controls_frame.grid(row=0, column=0, sticky='nsew')
		self.show_frame(ACSControlsFrame)
		
		self.frames[ACSPlotsFrame]=plots_frame
		plots_frame.grid(row=0, column=1, sticky='nsew')
		self.show_frame(ACSPlotsFrame)

		self.frames[ACSWavefrontFrame]=wavefront_frame
		wavefront_frame.grid(row=0, column=1, sticky='nsew')
		self.show_frame(ACSWavefrontFrame)
		
		self.frames[ACSImageFrame]=image_frame
		image_frame.grid(row=0, column=1, sticky='nsew')
		self.show_frame(ACSImageFrame)
		
		self.protocol("WM_DELETE_WINDOW", self.safe_quit)
		
		self.geometry('1250x750')

		self.filepath = ''

	
	def show_frame(self, cont):
		frame=self.frames[cont]
		frame.tkraise()


	def open_file(self):
		self.filepath = askopenfilename(title = 'Choose a file')

		if(len(self.filepath)==0):
			return None

		self.images_array, self.time_list, self.version, self.bitdepth = wfs.read_file(self.filepath)
		self.frames[ACSImageFrame].initialize_new_data(self.images_array[0], len(self.images_array), self.filepath)


	def get_parameters(self):
		control_frame = self.frames[ACSControlsFrame]
		self.path_length = float(control_frame.path_length_entry.get())
		self.aoi_size = int(control_frame.aoisize_entry.get())
		self.focal_length = float(control_frame.focallength_entry.get())*1e-3
		self.pixel_size = float(control_frame.pixelsize_entry.get())*1e-6
		self.wavelength = float(control_frame.wavelength_entry.get())*1e-9
		self.magnification = float(control_frame.magnification_entry.get())
		self.calculate_turbulence = bool(control_frame.calculate_turbulence.get())
		self.reconstruct = bool(control_frame.reconstruct_wavefront.get())


	def process_data(self):
		if(len(self.filepath)==0):
			return None
		image_frame = self.frames[ACSImageFrame]
		plots_frame = self.frames[ACSPlotsFrame]
		self.continue_processing = True

		self.get_parameters()
		file_processor = wfs.process_file(self.filepath, self.aoi_size, self.focal_length, self.pixel_size, self.wavelength, 
			self.magnification, calculate_turbulence=self.calculate_turbulence, reconstruct=self.reconstruct)

		if(not self.calculate_turbulence and not self.reconstruct):
			return None

		for frame in file_processor:
			status, return_values = frame
			if(status=='aoi locations'):
				aoi_locations = np.array([[int(val) for val in line.split(',')] for line in return_values])
				summed_array = np.sum(self.images_array, axis=0)
				image_frame.status_label.config(text='Status: Finding reference spots...')
				image_frame.add_aoi_locations(aoi_locations, self.aoi_size)
				image_frame.plot_data(summed_array)
			elif(status=='framenum'):
				image_frame.status_label.config(text='Status: Processing frames...')
				image_frame.update_framenum(return_values)
			else:
				if(status=='turbulence'):
					self.r0_full, self.r0_individual, self.ssf, self.ssfx, self.ssfy = return_values 
					self.frames[ACSPlotsFrame].set_slope_struct_functions(self.ssf, self.ssfx, self.ssfy, self.r0_full, self.r0_individual)
					self.frames[ACSPlotsFrame].plot_data()
				elif(status=='wavefront'):
					self.wavefronts_list = return_values
					self.frames[ACSWavefrontFrame].initialize_new_data(self.wavefronts_list)
					self.frames[ACSWavefrontFrame].plot_data()
				elif(status=='both'):
					self.r0_full, self.r0_individual, self.ssf, self.ssfx, self.ssfy, self.wavefronts_list = return_values
					self.frames[ACSPlotsFrame].set_slope_struct_functions(self.ssf, self.ssfx, self.ssfy, self.r0_full, self.r0_individual)
					self.frames[ACSPlotsFrame].plot_data()
					self.frames[ACSWavefrontFrame].initialize_new_data(self.wavefronts_list)
					self.frames[ACSWavefrontFrame].plot_data()
			
			
			if(not self.continue_processing):
				break
		image_frame.status_label.config(text='Status: Processing Complete')

	def thread_process_data(self):
		threading.Thread(target=self.process_data).start()


	def cancel_data_processing(self):
		self.continue_processing = False
		self.frames[ACSImageFrame].status_label.config(text='Status: Processing cancelled')




class ACSControlsFrame(tk.Frame):
	
	def __init__(self, parent, controller, imageframe, plotsframe):
		tk.Frame.__init__(self, parent)
		
		self.controller=controller
		self.image_frame=imageframe
		self.plots_frame=plotsframe
		self.grid_rowconfigure(0, weight=1)
		self.grid_columnconfigure(0, weight=1)
		self.current_image='Data'
		
		self.open_file_button=tk.Button(self, text='Open File', 
										command = self.controller.open_file, height=2, width=20)
		self.open_file_button.pack(side=tk.TOP, padx=30, pady=10)
		
		
		self.open_file_button=tk.Button(self, text='Process Data', 
										command = self.controller.thread_process_data, height=2, width=20)
		self.open_file_button.pack(side=tk.TOP, padx=30, pady=10)
		
		
		self.cancel_data_processing_button=tk.Button(self, text='Cancel Data Processing', 
										command= self.image_frame.cancel_data_processing, height=2, width=20)
		self.cancel_data_processing_button.pack(side=tk.TOP, padx=30, pady=10)
		
		
		self.toggle_image_button=tk.Button(self, text='Show Plots', 
										command= lambda: self.toggle_image(), height=2, width=20)
		self.toggle_image_button.pack(side=tk.TOP, padx=30, pady=10)

		self.calculate_turbulence = tk.BooleanVar()
		self.calculate_turbulence.set(False)
		self.calculate_turbulence_checkbutton=tk.Checkbutton(self, text='Calculate Turbulence', variable=self.calculate_turbulence)
		self.calculate_turbulence_checkbutton.pack(side=tk.TOP, padx=10, pady=(10,0))

		self.reconstruct_wavefront = tk.BooleanVar()
		self.reconstruct_wavefront.set(False)
		self.reconstruct_wavefront_checkbutton=tk.Checkbutton(self, text='Reconstruct Wavefront', variable=self.reconstruct_wavefront)
		self.reconstruct_wavefront_checkbutton.pack(side=tk.TOP, padx=10, pady=(10,0))

		
		self.path_length_label=tk.Label(self, text='Path Length (meters)')
		self.path_length_label.pack(side=tk.TOP, padx=10, pady=(10,0))
		self.path_length_entry=tk.Entry(self)
		self.path_length_entry.insert(0, '170')
		self.path_length_entry.pack(side=tk.TOP, padx=10, pady=0)
		
		
		self.wavelength_label=tk.Label(self, text='Wavelength (nm)')
		self.wavelength_label.pack(side=tk.TOP, padx=10, pady=(10,0))
		self.wavelength_entry=tk.Entry(self)
		self.wavelength_entry.insert(0, '640')
		self.wavelength_entry.pack(side=tk.TOP, padx=10, pady=0)


		self.magnification_label=tk.Label(self, text='Magnification')
		self.magnification_label.pack(side=tk.TOP, padx=10, pady=(10,0))
		self.magnification_entry=tk.Entry(self)
		self.magnification_entry.insert(0, '40')
		self.magnification_entry.pack(side=tk.TOP, padx=10, pady=0)


		self.aoisize_label=tk.Label(self, text='AOI Size (pixels)')
		self.aoisize_label.pack(side=tk.TOP, padx=10, pady=(10,0))
		self.aoisize_entry=tk.Entry(self)
		self.aoisize_entry.insert(0, '20')
		self.aoisize_entry.pack(side=tk.TOP, padx=10, pady=0)


		self.focallength_label=tk.Label(self, text='Lenslet Focal Length (mm)')
		self.focallength_label.pack(side=tk.TOP, padx=10, pady=(10,0))
		self.focallength_entry=tk.Entry(self)
		self.focallength_entry.insert(0, '6.7')
		self.focallength_entry.pack(side=tk.TOP, padx=10, pady=0)

		self.pixelsize_label=tk.Label(self, text='Pixel Size (microns)')
		self.pixelsize_label.pack(side=tk.TOP, padx=10, pady=(10,0))
		self.pixelsize_entry=tk.Entry(self)
		self.pixelsize_entry.insert(0, '7.4')
		self.pixelsize_entry.pack(side=tk.TOP, padx=10, pady=0)


		
		self.datetime_label=tk.Label(self, text='')
		self.datetime_label.pack(side=tk.TOP, padx=10, pady=10)
		
		self.version_label=tk.Label(self, text='')
		self.version_label.pack(side=tk.TOP, padx=10, pady=10)
		
		self.bitdepth_label=tk.Label(self, text='')
		self.bitdepth_label.pack(side=tk.TOP, padx=10, pady=10)


	def toggle_image(self):
		if(self.current_image=='Data'):
			self.current_image='Plots'
			self.toggle_image_button.config(text='Show Wavefront')
			self.controller.show_frame(ACSPlotsFrame)
		elif(self.current_image=='Plots'):
			self.current_image='Wavefront'
			self.toggle_image_button.config(text='Show Data')
			self.controller.show_frame(ACSWavefrontFrame)
		elif(self.current_image=='Wavefront'):
			self.current_image='Data'
			self.toggle_image_button.config(text='Show Plots')
			self.controller.show_frame(ACSImageFrame)



class ACSImageFrame(tk.Frame):
	
	def __init__(self, parent, controller):

		tk.Frame.__init__(self, parent)
		self.controller=controller
		
		self.control_frame=None
		self.grid_rowconfigure(0, weight=1)
		self.grid_columnconfigure(0, weight=1)
		
		self.filepath='No file selected yet'
		self.filepath_label=tk.Label(self, text='Filepath: '+self.filepath)
		self.filepath_label.pack(side=tk.TOP, padx=10, pady=10)
		
		self.current_framenum=0
		self.total_framenum=0
		self.framenum_label=tk.Label(self, text='Frame: '+str(self.current_framenum)+'/'+str(self.total_framenum))
		self.framenum_label.pack(side=tk.TOP, padx=10, pady=10)
		
		self.status_label=tk.Label(self, text='Status: Ready')
		self.status_label.pack(side=tk.TOP, padx=10, pady=10)
		
		self.f = Figure(figsize=(2,2), dpi=100)
		self.a = self.f.add_subplot(111)
		self.canvas = FigureCanvasTkAgg(self.f, self)
		self.canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
		self.canvas.draw()
		self.im=None

		self.aoi_locations = []
		self.aoi_size = 20
		
	
	def set_control_frame(self, control_frame):
		self.control_frame=control_frame


	def plot_data(self, image):
		self.a.clear()
		self.im = self.a.imshow(image)
		for aoi in self.aoi_locations:
			x, y = aoi[0], aoi[1]
			self.a.add_patch(patches.Circle((x+int(self.aoi_size/2),y+int(self.aoi_size/2)),radius=2,edgecolor='r',facecolor='r'))
			self.a.add_patch(patches.Rectangle((x,y), self.aoi_size, self.aoi_size, edgecolor='c', facecolor='none'))
		self.canvas.draw()


	def update_framenum(self, framenum):
		self.current_framenum = framenum
		self.framenum_label.config(text='Frame: '+str(self.current_framenum+1)+'/'+str(self.total_framenum))
		self.canvas.draw()


	def update_filepath(self, filepath):
		self.filepath = filepath
		self.filepath_label.config(text='Filepath: '+self.filepath)


	def initialize_new_data(self, image, numframes, filepath):
		self.plot_data(image)
		self.total_framenum = numframes
		self.update_filepath(filepath)
		self.update_framenum(self.current_framenum)


	def add_aoi_locations(self, aoi_locations, aoi_size):
		self.aoi_locations = aoi_locations
		self.aoi_size = aoi_size


	def cancel_data_processing(self):
		self.controller.cancel_data_processing()



class ACSPlotsFrame(tk.Frame):
	
	def __init__(self, parent, controller):
		tk.Frame.__init__(self, parent)
		
		self.controller=controller
		self.grid_rowconfigure(0, weight=1)
		self.grid_columnconfigure(0, weight=1)
		
		self.title='Slope Structure Functions'
		self.title_label=tk.Label(self, text=self.title)
		self.title_label.pack(side=tk.TOP, padx=10, pady=10)
		
		
		self.plot_toggle_button=tk.Button(self, text='Show XY SSF', 
										command=lambda: self.toggle_plots(), height=2, width=20)
		self.plot_toggle_button.pack(side=tk.TOP, padx=30, pady=10)
		
		
		self.f = Figure(figsize=(2,2), dpi=100)
		self.a = self.f.add_subplot(111)
		
		self.canvas = FigureCanvasTkAgg(self.f, self)
		self.canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
		self.canvas.draw()
		
		self.im = None

		self.current_plot = 'Full'
		self.ssf = []
		self.ssfx = []
		self.ssfy = []
		self.r0_full = None
		self.r0_individual = None
		
	
	
	def set_control_frame(self, controlframe):
		self.controlframe=controlframe
		return None
	

	def set_slope_struct_functions(self, ssf, ssfx, ssfy, r0_full, r0_individual):
		self.ssf = ssf
		self.ssfx = ssfx
		self.ssfy = ssfy
		self.r0_full = r0_full
		self.r0_individual = r0_individual
	

	def plot_data(self):
		self.a.clear()

		d_sub = self.controller.pixel_size*self.controller.aoi_size*self.controller.magnification

		def functional_form_full(r, r0):
			val = np.multiply(6.88/(d_sub**2)*(d_sub/r0)**(5/3)*1.32, np.add(np.divide(np.square(r),np.add(5.4,np.square(r))),np.multiply(0.557, np.power(r, 0.213))))
			return val

		def functional_form_individual(r, r0):
			val = np.multiply(6.88/(d_sub**2)*(d_sub/r0)**(5/3)*0.77, np.add(np.divide(np.square(r),np.add(1,np.square(r))),np.multiply(0.438, np.power(r, 0.2555))))
			return val        
		if(self.current_plot=='Full'):
			r0 = self.r0_full
			functional_form = functional_form_full
			ssf = self.ssf
			r = np.linspace(0,ssf[-1,0],1000)
			self.a.plot(ssf[:,0], ssf[:,1], 'r--', label='Data')
			self.a.plot(r, functional_form(r, r0), 'r-', label='Fit')
			self.a.legend(prop={'size':15})
			self.a.set_xlabel(r'$\Delta r$')
			self.a.set_ylabel(r'$D_s$')
			self.canvas.draw()
		elif(self.current_plot=='XY'):
			r0 = self.r0_individual
			functional_form = functional_form_individual
			ssfx, ssfy = self.ssfx, self.ssfy
			r = np.linspace(0,ssfx[-1,0],1000)
			self.a.plot(ssfx[:,0], ssfx[:,1], 'r--', label='X Data')
			self.a.plot(ssfy[:,0], ssfy[:,1], 'g--', label='Y Data')
			self.a.plot(r, functional_form(r, r0), 'b-', label='Fit')
			self.a.legend(prop={'size':15})
			self.a.set_xlabel(r'$\Delta r$')
			self.a.set_ylabel(r'$D_{x,y}$')
			self.canvas.draw()

	def toggle_plots(self):
		if(self.current_plot=='Full'):
			self.current_plot='XY'
			self.plot_toggle_button.config(text='Show Full SSF')
		elif(self.current_plot=='XY'):
			self.current_plot='Full'
			self.plot_toggle_button.config(text='Show XY SSF')
		self.plot_data()



class ACSWavefrontFrame(tk.Frame):
	
	def __init__(self, parent, controller):

		tk.Frame.__init__(self, parent)
		self.controller=controller
		
		self.control_frame=None
		self.grid_rowconfigure(0, weight=1)
		self.grid_columnconfigure(0, weight=1)
		
		self.current_framenum=0
		self.total_framenum=0
		self.framenum_label=tk.Label(self, text='Frame: '+str(self.current_framenum)+'/'+str(self.total_framenum))
		self.framenum_label.pack(side=tk.TOP, padx=10, pady=10)

		
		self.f = Figure(figsize=(2,2), dpi=100)
		self.a = self.f.add_subplot(111)
		self.canvas = FigureCanvasTkAgg(self.f, self)
		self.canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
		self.canvas.draw()
		self.im=None

		self.wavefronts_list = []

		self.whitespace_label=tk.Label(self, text='                                       ')
		self.whitespace_label.pack(side=tk.RIGHT, padx=30, pady=10)

		self.next_frame_button=tk.Button(self, text='Next Frame', 
								command= self.plot_next, height=2, width=20)
		self.next_frame_button.pack(side=tk.RIGHT, padx=30, pady=10)


		self.whitespace_label=tk.Label(self, text='                                       ')
		self.whitespace_label.pack(side=tk.LEFT, padx=30, pady=10)


		self.last_frame_button=tk.Button(self, text='Last Frame', 
								command= self.plot_last, height=2, width=20)
		self.last_frame_button.pack(side=tk.LEFT, padx=30, pady=10)
		
	
	def set_control_frame(self, control_frame):
		self.control_frame=control_frame


	def plot_data(self):
		self.a.clear()
		self.im = self.a.imshow(self.wavefronts_list[self.current_framenum])
		self.canvas.draw()


	def update_framenum(self, framenum):
		self.current_framenum = framenum
		self.framenum_label.config(text='Frame: '+str(self.current_framenum+1)+'/'+str(self.total_framenum))
		self.canvas.draw()


	def initialize_new_data(self, wavefronts_list):
		self.wavefronts_list = wavefronts_list
		if(len(self.wavefronts_list)==0):
			return None
		self.current_framenum = 0
		self.plot_data()
		self.total_framenum = len(self.wavefronts_list)
		self.update_framenum(self.current_framenum)


	def plot_next(self):
		if(self.current_framenum<self.total_framenum-1):
			self.current_framenum+=1
			self.plot_data()
			self.update_framenum(self.current_framenum)


	def plot_last(self):
		if(self.current_framenum>0):
			self.current_framenum-=1
			self.plot_data()
			self.update_framenum(self.current_framenum)




if __name__=='__main__':
	app = ACSDataApp()
	app.mainloop()