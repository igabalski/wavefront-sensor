

# In[1]:


# Example for loading .dat binary files.
#
# Date and time stamps come from the .dat files names and corresponds to
# the date and time at which the image acquisition begain.
#
# Follow the logic in the script to see what the order of the data saved in
# the files and each parameters data type.
#
# This script must be in the same directory as the .dat file(s) to execute.
#
# 13-Feb-2018, D. Whitley, Script creation
# 1-May-2018, I. Gabalski, script modification


import numpy as np
from scipy.ndimage.measurements import center_of_mass
from scipy.optimize import curve_fit
import random
import tkinter as tk
from tkinter.filedialog import askopenfilename, asksaveasfilename
import os
import warnings
warnings.filterwarnings("ignore")
import time
from itertools import combinations
import math

from tkinter import ttk
import matplotlib
import matplotlib.patches as patches
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib.backend_bases import key_press_handler
from mpl_toolkits.axes_grid1 import make_axes_locatable

import codecs
import threading
import struct


# In[2]:


class ACSDataApp(tk.Tk):
    
    def __init__(self, *args, **kwargs):
        
        tk.Tk.__init__(self, *args, **kwargs)
        
        container=tk.Frame(self)
        container.pack(side='top', fill='both', expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)
        container.grid_columnconfigure(1, weight=15)
        
        self.wm_title('ACS Data Analysis')
        
        self.frames={}
        clickerframe=ACSDataFrame(container, self)
        plotsframe=ACSPlotsFrame(container, self)
        ssplotsframe=ACSSlopeStructPlotFrame(container, self)
        controlframe=ACSDataControls(container, self, clickerframe, plotsframe)
        clickerframe.set_control_frame(controlframe)
        plotsframe.set_control_frame(controlframe)
        ssplotsframe.set_control_frame(controlframe)
        
        
        self.frames[ACSDataControls]=controlframe
        controlframe.grid(row=0, column=0, sticky='nsew')
        self.show_frame(ACSDataControls)
        
        self.frames[ACSPlotsFrame]=plotsframe
        plotsframe.grid(row=0, column=1, sticky='nsew')
        self.show_frame(ACSPlotsFrame)
        
        self.frames[ACSSlopeStructPlotFrame]=ssplotsframe
        ssplotsframe.grid(row=0, column=1, sticky='nsew')
        self.show_frame(ACSSlopeStructPlotFrame)
        
        self.frames[ACSDataFrame]=clickerframe
        clickerframe.grid(row=0, column=1, sticky='nsew')
        self.show_frame(ACSDataFrame)
        
        
        
        self.geometry('1250x750')
        
        
        
    
    def show_frame(self, cont):
        frame=self.frames[cont]
        frame.tkraise()
    
    def safe_quit(self):
        self.frames[ACSDataFrame].cancel_data_processing()
        time.sleep(0.3)
        self.quit()
        self.destroy()


# In[3]:


class ACSDataControls(tk.Frame):
    
    def __init__(self, parent, controller, imageframe, plotsframe):
        tk.Frame.__init__(self, parent)
        
        self.controller=controller
        self.imageframe=imageframe
        self.plotsframe=plotsframe
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        
        self.open_file_button=tk.Button(self, text='Open and Process Data', 
                                        command= self.imageframe.thread_open_file, height=2, width=20)
        self.open_file_button.pack(side=tk.TOP, padx=30, pady=10)
        
        
#         self.references_button=tk.Button(self, text='Calculate references', 
#                                         command= self.imageframe.calculate_reference_spots, height=2, width=20)
#         self.references_button.pack(side=tk.TOP, padx=30, pady=10)
        
        
#         self.centroids_button=tk.Button(self, text='Calculate centroids', 
#                                         command= self.imageframe.calculate_centroids, height=2, width=20)
#         self.centroids_button.pack(side=tk.TOP, padx=30, pady=10)
        
        
#         self.calculate_cn2_button=tk.Button(self, text='Calculate Cn^2', 
#                                         command= self.imageframe.calculate_cn2, height=2, width=20)
#         self.calculate_cn2_button.pack(side=tk.TOP, padx=30, pady=10)
        
#         self.generate_plots_button=tk.Button(self, text='Generate Plots', 
#                                         command= self.plotsframe.plot_data, height=2, width=20)
#         self.generate_plots_button.pack(side=tk.TOP, padx=30, pady=10)
        
        
        self.cancel_data_processing_button=tk.Button(self, text='Cancel Data Processing', 
                                        command= self.imageframe.cancel_data_processing, height=2, width=20)
        self.cancel_data_processing_button.pack(side=tk.TOP, padx=30, pady=10)
        
        
        self.show_data_button=tk.Button(self, text='Show Data', 
                                        command= lambda: self.controller.show_frame(ACSDataFrame), height=2, width=20)
        self.show_data_button.pack(side=tk.TOP, padx=30, pady=10)
        
        self.show_plots_button=tk.Button(self, text='Show Cn2 Plots', 
                                        command= lambda: self.controller.show_frame(ACSPlotsFrame), height=2, width=20)
        self.show_plots_button.pack(side=tk.TOP, padx=30, pady=10)
        
        self.show_ssplots_button=tk.Button(self, text='Show Slope Structure Plots', 
                                        command= lambda: self.controller.show_frame(ACSSlopeStructPlotFrame), height=2, width=20)
        self.show_ssplots_button.pack(side=tk.TOP, padx=30, pady=10)
        
        self.threshold_label=tk.Label(self, text='Ref. Threshold (% max)')
        self.threshold_label.pack(side=tk.TOP, padx=10, pady=10)
        self.threshold_entry=tk.Entry(self)
        self.threshold_entry.insert(0, '30')
        self.threshold_entry.pack(side=tk.TOP, padx=10, pady=0)
        
        
        self.aoi_label=tk.Label(self, text='AOI Threshold')
        self.aoi_label.pack(side=tk.TOP, padx=10, pady=10)
        self.aoi_entry=tk.Entry(self)
        self.aoi_entry.insert(0, '4')
        self.aoi_entry.pack(side=tk.TOP, padx=10, pady=0)
        
        
        self.path_length_label=tk.Label(self, text='Path Length (meters)')
        self.path_length_label.pack(side=tk.TOP, padx=10, pady=10)
        self.path_length_entry=tk.Entry(self)
        self.path_length_entry.insert(0, '1000')
        self.path_length_entry.pack(side=tk.TOP, padx=10, pady=0)
        
        
        self.wavelength_label=tk.Label(self, text='Wavelength (nm)')
        self.wavelength_label.pack(side=tk.TOP, padx=10, pady=10)
        self.wavelength_entry=tk.Entry(self)
        self.wavelength_entry.insert(0, '640')
        self.wavelength_entry.pack(side=tk.TOP, padx=10, pady=0)
        
        
        self.datetime_label=tk.Label(self, text='')
        self.datetime_label.pack(side=tk.TOP, padx=10, pady=10)
        
        self.version_label=tk.Label(self, text='')
        self.version_label.pack(side=tk.TOP, padx=10, pady=10)
        
        self.bitdepth_label=tk.Label(self, text='')
        self.bitdepth_label.pack(side=tk.TOP, padx=10, pady=10)


# In[7]:


class ACSDataFrame(tk.Frame):
    
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller=controller
        self.timeIn = []
        self.dataArray = []
        self.imagesArray = []
        
        self.controlframe=None
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        
        self.current_filepath='No file selected yet'
        self.filepath_label=tk.Label(self, text='Filepath: '+self.current_filepath)
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
        
        self.imgFrameData=None
        self.im=None
        self.animating=False
        
        cid=self.canvas.mpl_connect('button_press_event', self.on_key_event)
        
        self.summed_array=[]
        self.summed_array_unthresholded=[]
        self.reference_spots=[]
        self.sorted_reference_spots=[]
        self.centroids_list=[]
        self.displacements_list=[]
        self.cn2_calculated=[]
        self.r0_calculated=[]
        self.slopediffs_x=[]
        self.slopediffs_y=[]
        self.r0s=[]
        self.num_framesprocessed=0
        self.num_framesread=0
        
        
        self.process_data=False
        self.read_data=False
        self.magnification=2/0.05
        self.aperture_size=0.203
        self.pixel_size=7.4e-6
        self.subaperture_size=20
    
    
    def set_control_frame(self, controlframe):
        self.controlframe=controlframe
        return None
            
            
    
    #This method opens the file and stores the data in an array
    #Look here to find file format
    #need to put images into several arrays in order to avoid issues with large files
    def thread_open_file(self):
        threading.Thread(target=self.open_file).start()
    
    
    
    
    def open_file(self):
        self.status_label.config(text='Status: Reading File...')
        self.timeIn = []
        self.dataArray = []
        self.imagesArray = []
        self.cn2_calculated=[]
        self.r0_calculated=[]
        self.slopediffs_x=[]
        self.slopediffs_y=[]
        self.r0s=[]
        self.num_framesprocessed=0
        self.num_framesread=0
        try:
            filepath=askopenfilename(title="Choose a file", 
                                     initialdir='F:/Black Dart 2017/ACS Data/05_19_2017/data2017_05_19_00_01_01')
            self.filepath_label.config(text=filepath)
            fileName=filepath.split('/')[-1]

            year = fileName[4:8]
            month = fileName[9:11]
            day = fileName[12:14]
            hour = fileName[15:17]
            minute = fileName[18:20]
            second = fileName[21:23]
            
        except Exception as e:
            print(e)
        
        try:
            if len(filepath)>0:
                self.controlframe.datetime_label.config(text="%s/%s/%s %s:%s:%s UTC" %(month,day,year,hour,minute,second))
            self.read_data=True
            with open(filepath,"rb") as f:

                version = np.fromfile(f,count=1,dtype=np.double)
                bitdepth = np.fromfile(f,count=1,dtype=np.int32)
                self.controlframe.version_label.config(text="version = %s" % (version))
                self.controlframe.bitdepth_label.config(text="bitdepth = %s" % (bitdepth))

                self.process_data=True
                
                for i in range(500):
                    frameID = np.fromfile(f,count=1,dtype=np.uint64)
                    if(len(frameID)==0):
                        self.read_data=False
                        break
                    thisTime = np.fromfile(f,count=1,dtype=np.double)
                    self.timeIn.append(thisTime)
                    width = int(np.fromfile(f,count=1,dtype=np.uint32))
                    height = int(np.fromfile(f,count=1,dtype=np.uint32))
                    bufSize = np.fromfile(f,count=1,dtype=np.uint32) 
                    #buffer size is size of a single frame in bytes. Frame reads out as a stream of pixel values.
                    
                    thisArray = np.fromfile(f,count=int(width*height),dtype=np.uint16)
                    if(len(thisArray) == (width*height)):
                        self.dataArray.append(thisArray)
                        thisImageArray = thisArray.reshape([height,width])
                        self.imagesArray.append(thisImageArray)
                        self.num_framesread+=1
                    else:
                        print("Frame %s Skipped!!!" % (frameID))
                self.a.clear()
                self.im=self.a.imshow(self.imagesArray[0])
                self.canvas.draw()
                
                self.framenum_label.config(text='Frame: '+str(1)+'/'+str(len(self.imagesArray)))
                self.calculate_reference_spots()
                self.process_data_method()
                
                
                while(self.read_data): # this loops through the file until it hits EOF
                    self.imagesArray=[]
                    num_framestoread=1000
                    for i in range(num_framestoread):
                        frameID = np.fromfile(f,count=1,dtype=np.uint64)
                        if(len(frameID)==0):
                            self.read_data=False
                            break  
                        thisTime = np.fromfile(f,count=1,dtype=np.double)
                        self.timeIn.append(thisTime)
                        width = int(np.fromfile(f,count=1,dtype=np.uint32))
                        height = int(np.fromfile(f,count=1,dtype=np.uint32))
                        bufSize = np.fromfile(f,count=1,dtype=np.uint32) 
                        #buffer size is size of a single frame in bytes. Frame reads out as a stream of pixel values.
                        
                        thisArray = np.fromfile(f,count=int(width*height),dtype=np.uint16)
                        if(len(thisArray) == (width*height)):
                            self.dataArray.append(thisArray)
                            thisImageArray = thisArray.reshape([height,width])
                            self.imagesArray.append(thisImageArray)
                            self.num_framesread+=1
                        else:
                            print("Frame %s Skipped!!!" % (frameID))
                    if(not self.read_data):
                        break
                    else:
                        self.process_data_method()
                    self.controller.frames[ACSSlopeStructPlotFrame].r0_framenumber+=num_framestoread
                
        except Exception as e:
            print(e)
        print('Done Processing')
        print(np.mean(self.cn2_calculated), np.std(self.cn2_calculated))
    
    
    def on_key_event(self, event):
        if event.inaxes is not None:
            if event.button==1:
                key_press_handler(event, self.canvas)
                if(self.current_framenum<self.num_framesread):
                    self.im.set_data(self.imagesArray[self.current_framenum])
                    self.canvas.draw()
                    self.framenum_label.config(text='Frame: '+str(self.current_framenum)+'/'+str(self.total_framenum))
    
    
    
    
    
    def calculate_reference_spots(self):
        self.status_label.config(text='Status: Calculating Reference Spots...')
        self.reference_spots=[]
        #sum and threshold all frames to average
        threshold=int(self.controlframe.threshold_entry.get())/100
        self.summed_array_unthresholded=np.sum(self.imagesArray, axis=0)
        self.summed_array=np.array([i for i in self.summed_array_unthresholded])
        maxval=np.amax(self.summed_array)
        pixel_threshold=maxval*threshold
        self.summed_array[self.summed_array<pixel_threshold]=0
        
        
        self.im.set_data(self.summed_array)
        self.canvas.draw()
        test_window_size=15
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
        
        local_maxima[:,0:int(test_window_size/2)]=0
        local_maxima[:,xlim-int(test_window_size/2):xlim]=0
        local_maxima[0:int(test_window_size/2),:]=0
        local_maxima[ylim-int(test_window_size/2):ylim,:]=0
        
        
        #add found reference spots to array, and red circles to image for testing
        window_size=17
        window_offset=8 #(window_size-1)/2
        for y in range(len(local_maxima)):
            for x in range(len(local_maxima[y])):
                if(local_maxima[y][x]==1):
                    self.a.add_patch(patches.Circle((x,y),radius=2,edgecolor='r',facecolor='r'))
                    self.a.add_patch(patches.Rectangle((x-9,y-9),17,17, edgecolor='c', facecolor='none'))
                    
                    rel_cent_x, rel_cent_y=center_of_mass(self.summed_array_unthresholded[y-window_offset:y-window_offset+window_size,x-window_offset:x-window_offset+window_size])
                    cent_x=rel_cent_x+x-window_offset
                    cent_y=rel_cent_y+y-window_offset
                    
                    self.reference_spots.append([x,y, cent_x, cent_y])
                if(not self.process_data):
                    break
            if(not self.process_data):
                break
        self.canvas.draw()
        
    
    
    def calculate_centroids(self):
        self.displacements_list=[]
        window_size=17
        window_offset=8 #(window_size-1)/2
        framenum=0
        dropped_aoi_threshold_factor=int(self.controlframe.aoi_entry.get()) #for excluding AOIs for which maxval<meanval*this factor
        
        for frame in self.imagesArray:
            centroids=[]
            displacements=[]
            
            for spot in self.reference_spots:
                spot_x=spot[0]
                spot_y=spot[1]
                spot_cent_x=spot[2]
                spot_cent_y=spot[3]
                aoi=np.array(frame[spot_y-window_offset:spot_y-window_offset+window_size,spot_x-window_offset:spot_x-window_offset+window_size])
                maxval=np.amax(aoi)
                threshold=maxval*0.4
                aoi[aoi<threshold]=0
                rel_cent_x, rel_cent_y=center_of_mass(aoi)
                cent_x=rel_cent_x+spot_cent_x-window_offset
                cent_y=rel_cent_y+spot_cent_y-window_offset
                
                maxval=np.amax(frame[spot_y-window_offset:spot_y-window_offset+window_size,spot_x-window_offset:spot_x-window_offset+window_size])
                meanval=frame[spot_y-window_offset:spot_y-window_offset+window_size,spot_x-window_offset:spot_x-window_offset+window_size].mean()
                if((not math.isnan(rel_cent_x)) and (maxval>meanval*dropped_aoi_threshold_factor)):
                    #this appends ROI locations in subaperture-normalized units, displacements in meters at sensor
                    displacements.append([spot_x/self.subaperture_size, spot_y/self.subaperture_size,(spot_x-cent_x)*self.pixel_size, (spot_y-cent_y)*self.pixel_size])
                if(not self.process_data):
                    break
            if(not self.process_data):
                break
            
            self.displacements_list.append(displacements)
            framenum+=1
            
            xdiff_mean=np.nanmean(np.array(displacements)[:,2])
            ydiff_mean=np.nanmean(np.array(displacements)[:,3])
            self.slopediffs_x.append(xdiff_mean)
            self.slopediffs_y.append(ydiff_mean)
            
            
            self.framenum_label.config(text='Frame: '+str(framenum)+'/'+str(len(self.imagesArray)))
            

        
    
    def calculate_cn2(self):
        self.slope_struct_functions_x=[]
        self.slope_struct_functions_y=[]
        rbin_size_in_pixels=20
        rbin_size=rbin_size_in_pixels/self.subaperture_size
        num_rbins=int(np.sqrt(len(self.imagesArray[0])**2+len(self.imagesArray[0][0])**2)/rbin_size_in_pixels)
        framenum=0
        
        
        wavenumber=2*np.pi/(int(self.controlframe.wavelength_entry.get())*10**(-9))
        focal_length=0.0067
        
        framenum=0
        #frame is of format (x, y, dx, dy) all values normalized to subaperture size
        for frame in self.displacements_list:
            self.framenum_label.config(text='Frame: '+str(self.num_framesprocessed)+'/'+str(self.num_framesread))
            
            #create list of pairs of AOIs
            frame_pairs=np.array(list(combinations(frame, 2)))
            
            #(k/(f*m))^2 factor that pulls out of slope structure function
            factor=(wavenumber/focal_length/self.magnification)**2
            
            unbinned_radii=np.sqrt(np.add(np.square(np.subtract(frame_pairs[:,0,0], frame_pairs[:,1,0])),np.square(np.subtract(frame_pairs[:,0,1], frame_pairs[:,1,1]))))
            unbinned_xslopes_squared=np.multiply(factor, np.square(np.subtract(frame_pairs[:,0,2],frame_pairs[:,1,2])))
            unbinned_yslopes_squared=np.multiply(factor, np.square(np.subtract(frame_pairs[:,0,3],frame_pairs[:,1,3])))
            slope_struct_func_x=np.array([np.nanmean(unbinned_xslopes_squared[np.logical_and(rbin<unbinned_radii,unbinned_radii<=rbin+rbin_size)]) for rbin in np.arange(0,num_rbins*rbin_size,rbin_size)])
            slope_struct_func_y=np.array([np.nanmean(unbinned_yslopes_squared[np.logical_and(rbin<unbinned_radii,unbinned_radii<=rbin+rbin_size)]) for rbin in np.arange(0,num_rbins*rbin_size,rbin_size)])
            
            
            slope_struct_func_x[np.isnan(slope_struct_func_x)]=0
            slope_struct_func_y[np.isnan(slope_struct_func_y)]=0
            
            slope_struct_func_x=np.trim_zeros(slope_struct_func_x, trim='b')
            slope_struct_func_y=np.trim_zeros(slope_struct_func_y, trim='b')
            
            slope_struct_func_x=slope_struct_func_x[:-1]
            slope_struct_func_y=slope_struct_func_y[:-1]
            
            fudge_factor=1
            slope_struct_func_x=np.multiply(fudge_factor,slope_struct_func_x)
            slope_struct_func_y=np.multiply(fudge_factor,slope_struct_func_y)
            
            self.slope_struct_functions_x.append(slope_struct_func_x)
            self.slope_struct_functions_y.append(slope_struct_func_y)
            if(not self.process_data):
                break
            
            self.status_label.config(text='Status: Calculating Slope Structure Functions...')
            self.a.clear()
            self.im=self.a.imshow(self.imagesArray[framenum])
            self.canvas.draw()
            self.num_framesprocessed+=1
            framenum+=1
            
            
        self.status_label.config(text='Status: Fitting Slope Structure Functions...')
        
        for i in range(len(self.slope_struct_functions_x)):
            self.framenum_label.config(text='Frame: '+str(i)+'/'+str(len(self.slope_struct_functions_x)))
            d_sub=self.pixel_size*self.subaperture_size*self.magnification
            def functional_form(r, r0):
                return 0.77*6.88/(d_sub**2)*(d_sub/r0)**(5/3)*(r**2/(1+r**2)+0.438*r**(0.2555))
            
            r_vals_x=[rbin_size*i+rbin_size/2 for i in range(len(self.slope_struct_functions_x[i]))]
            r_vals_y=[rbin_size*i+rbin_size/2 for i in range(len(self.slope_struct_functions_y[i]))]
            xparams, xcovariances=curve_fit(functional_form, r_vals_x, self.slope_struct_functions_x[i])
            yparams, ycovariances=curve_fit(functional_form, r_vals_y, self.slope_struct_functions_y[i])
            r0_x=xparams[0]
            r0_y=yparams[0]
            self.r0s.append((r0_x,r0_y))
            r_0=np.mean([r0_x,r0_y])
            self.r0_calculated.append(r_0)
            #calculate Cn^2 from r_0
            path_length=int(self.controlframe.path_length_entry.get())
            wavelength=int(self.controlframe.wavelength_entry.get())*10**(-9)
            cn2=wavelength**2/(4*np.pi**2)*(r_0/3)**(-5/3)/path_length
            self.cn2_calculated.append(cn2)
            if(not self.process_data):
                break
            
        
            
    
    def thread_process_data(self):
        threading.Thread(target=self.process_data_method).start()
    
    def process_data_method(self):
        self.process_data=True
        if(self.process_data):
            self.status_label.config(text='Status: Calculating centroids...')
            self.calculate_centroids()
        if(self.process_data):
            self.calculate_cn2()
        if(self.process_data):
            self.status_label.config(text='Status: Plotting data...')
            self.controller.frames[ACSPlotsFrame].plot_data()
            self.controller.frames[ACSSlopeStructPlotFrame].set_slope_struct_functions(self.slope_struct_functions_x, self.slope_struct_functions_y, self.r0_calculated, self.r0s)
            self.controller.frames[ACSSlopeStructPlotFrame].plot_data()
        if(self.process_data):
            self.status_label.config(text='Status: Processing Complete')
        else:
            self.status_label.config(text='Status: Processing Cancelled')
        self.process_data=False
    
    
    def cancel_data_processing(self):
        self.process_data=False
        self.read_data=False
        
        
        


# In[8]:


class ACSPlotsFrame(tk.Frame):
    
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        
        self.controller=controller
        self.controlframe=None
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        
        self.title='Cn2 Plots'
        self.title_label=tk.Label(self, text=self.title)
        self.title_label.pack(side=tk.TOP, padx=10, pady=10)
        
        self.last_frame_button=tk.Button(self, text='Plot Cn2', 
                                        command= lambda: self.plot_cn2(), height=2, width=20)
        self.last_frame_button.pack(side=tk.TOP, padx=20, pady=10)
        
        self.last_frame_button=tk.Button(self, text='Plot r0', 
                                        command= lambda: self.plot_r0(), height=2, width=20)
        self.last_frame_button.pack(side=tk.TOP, padx=20, pady=10)
        
        self.last_frame_button=tk.Button(self, text='Plot Avg Displacements', 
                                        command= lambda: self.plot_avg_displacements(), height=2, width=20)
        self.last_frame_button.pack(side=tk.TOP, padx=20, pady=10)
        
        self.whichplot='cn2'
        
        self.f = Figure(figsize=(2,2), dpi=100)
        self.a = self.f.add_subplot(111)
        
        self.canvas = FigureCanvasTkAgg(self.f, self)
        self.canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        self.canvas.draw()
        
        self.im=None
        
        self.slope_struct_functions_x=[]
        self.slope_struct_functions_y=[]
    
    
    def set_control_frame(self, controlframe):
        self.controlframe=controlframe
        return None
            
    def set_slope_struct_functions(self, slope_struct_func_x, slope_struct_func_y):
        self.slope_struct_functions_x=slope_struct_func_x
        self.slope_struct_functions_y=slope_struct_func_y
    
    def plot_data(self):
        if(self.whichplot=='cn2'):
            self.a.clear()
            self.x=self.a.semilogy(self.controller.frames[ACSDataFrame].cn2_calculated, label=r'$C_n^2$')
            self.a.legend(prop={'size':15})
            self.a.set_xlabel(r'Frame#')
            self.a.set_ylabel(r'$C_n^2 (m^{-\frac{2}{3}})$')
            cn2_min=np.amin(self.controller.frames[ACSDataFrame].cn2_calculated)
            cn2_max=np.amax(self.controller.frames[ACSDataFrame].cn2_calculated)
            #self.a.set_yticks(np.linspace(cn2_min,cn2_max,10))
            self.canvas.draw()
        elif(self.whichplot=='r0'):
            self.a.clear()
            self.x=self.a.plot(self.controller.frames[ACSDataFrame].r0_calculated, label=r'$r_0$')
            self.a.legend(prop={'size':15})
            self.a.set_xlabel(r'Frame#')
            self.a.set_ylabel(r'$r_0$')
            r0_min=np.amin(self.controller.frames[ACSDataFrame].r0_calculated)
            r0_max=np.amax(self.controller.frames[ACSDataFrame].r0_calculated)
            #self.a.set_yticks(np.linspace(r0_min,r0_max,10))
            self.canvas.draw()
        elif(self.whichplot=='avgd'):
            self.a.clear()
            self.x=self.a.plot(self.controller.frames[ACSDataFrame].slopediffs_x, label=r'$x$')
            self.x=self.a.plot(self.controller.frames[ACSDataFrame].slopediffs_y, label=r'$y$')
            zeros=[0 for i in range(len(self.controller.frames[ACSDataFrame].slopediffs_x))]
            self.x=self.a.plot(zeros, label=r'$0$')
            self.a.legend(prop={'size':15})
            self.a.set_xlabel(r'Frame#')
            self.a.set_ylabel(r'$Avg. Displacement$')
            self.canvas.draw()
    
    def plot_r0(self):
        self.whichplot='r0'
        self.plot_data()
    
    def plot_cn2(self):
        self.whichplot='cn2'
        self.plot_data()
        
    def plot_avg_displacements(self):
        self.whichplot='avgd'
        self.plot_data()
    

class ACSSlopeStructPlotFrame(tk.Frame):
    
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        
        self.controller=controller
        self.controlframe=None
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        
        self.title='Slope Structure Functions'
        self.title_label=tk.Label(self, text=self.title)
        self.title_label.pack(side=tk.TOP, padx=10, pady=10)
        
        
        self.last_frame_button=tk.Button(self, text='Last Frame', 
                                        command= lambda: self.show_last_frame(), height=2, width=20)
        self.last_frame_button.pack(side=tk.TOP, padx=30, pady=10)
        
        
        self.next_frame_button=tk.Button(self, text='Next Frame', 
                                        command= lambda: self.show_next_frame(), height=2, width=20)
        self.next_frame_button.pack(side=tk.TOP, padx=30, pady=10)
        
        
        self.f = Figure(figsize=(2,2), dpi=100)
        self.a = self.f.add_subplot(111)
        
        self.canvas = FigureCanvasTkAgg(self.f, self)
        self.canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        self.canvas.draw()
        
        self.im=None
        
        self.slope_struct_functions_x=[]
        self.slope_struct_functions_y=[]
        self.r0=[]
        self.r0s=[]
        self.framenum=0
        self.r0_framenum=0
    
    
    def set_control_frame(self, controlframe):
        self.controlframe=controlframe
        return None
            
    def set_slope_struct_functions(self, slope_struct_funcs_x, slope_struct_funcs_y, r0_vals, r0s_vals):
        self.slope_struct_functions_x=slope_struct_funcs_x
        self.slope_struct_functions_y=slope_struct_funcs_y
        self.r0=r0_vals
        self.r0s=r0s_vals
    
    def plot_data(self):
        self.a.clear()
        rx=np.arange(len(self.slope_struct_functions_x[self.framenum]))
        ry=np.arange(len(self.slope_struct_functions_y[self.framenum]))
        r_form=np.linspace(0,35,10000)
        D_form_x=[0.77*6.88/((7.4e-6*20*40)**2)*(7.4e-6*20*40/(self.r0s[self.r0_framenum][0]))**(5/3)*(r**2/(1+r**2)+0.438*r**(0.2555)) for r in r_form]
        D_form_y=[0.77*6.88/((7.4e-6*20*40)**2)*(7.4e-6*20*40/(self.r0s[self.r0_framenum][1]))**(5/3)*(r**2/(1+r**2)+0.438*r**(0.2555)) for r in r_form]
        self.x=self.a.plot(rx, self.slope_struct_functions_x[self.framenum], label=r'$D_x$')
        self.x=self.a.plot(ry, self.slope_struct_functions_y[self.framenum], label=r'$D_y$')
        self.x=self.a.plot(r_form, D_form_x, label=r'$D_{formx}$')
        self.x=self.a.plot(r_form, D_form_y, label=r'$D_{formy}$')
        self.a.legend(prop={'size':15})
        self.a.set_xlabel(r'$\Delta r$')
        self.a.set_ylabel(r'$D_s$')
        self.canvas.draw()
        
        
    def show_last_frame(self):
        if(self.framenum>0):
            self.framenum-=1
            self.r0_framenum-=1
            self.plot_data()
            
    
    def show_next_frame(self):
        if(self.framenum<len(self.slope_struct_functions_x)-1):
            self.framenum+=1
            self.r0_framenum+=1
            self.plot_data()
            




    
        


# In[9]:

print('Notes:')
print('1. Slow step is filling unbinned_differences with values')
print('2. Single frame method not sensitive to tilt - spots move together implies difference is 0')
print('3. Same frame being saved several times - talk to Brett about Labview software')
app=ACSDataApp()
app.protocol("WM_DELETE_WINDOW", app.safe_quit)
app.mainloop()