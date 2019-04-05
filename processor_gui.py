import tkinter as tk
from tkinter.filedialog import askopenfilename, asksaveasfilename, askdirectory
import os
import time

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
        
        
        self.frames[ACSControlsFrame]=controls_frame
        controls_frame.grid(row=0, column=0, sticky='nsew')
        self.show_frame(ACSControlsFrame)
        
        self.frames[ACSPlotsFrame]=plots_frame
        plots_frame.grid(row=0, column=1, sticky='nsew')
        self.show_frame(ACSPlotsFrame)
        
        self.frames[ACSImageFrame]=image_frame
        image_frame.grid(row=0, column=1, sticky='nsew')
        self.show_frame(ACSImageFrame)
        
        self.protocol("WM_DELETE_WINDOW", self.safe_quit)
        
        self.geometry('1250x750')

    
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


    def process_data(self):
        image_frame = self.frames[ACSImageFrame]
        plots_frame = self.frames[ACSPlotsFrame]

        self.get_parameters()
        file_processor = wfs.process_file(self.filepath, self.aoi_size, self.focal_length, self.pixel_size, self.wavelength, self.magnification)
        for frame in file_processor:
            return_values = frame
            if(isinstance(return_values, list)):
                aoi_locations = np.array([[int(val) for val in line.split(',')] for line in return_values])
                ...update gui with aoi locations...
            elif(isinstance(return_values, tuple)):
                ...update gui with 'fitting ssf data' message...
            else:
                ...update gui with frame number...
        r0_full, r0_individual, ssf, ssfx, ssfy = return_values 




class ACSControlsFrame(tk.Frame):
    
    def __init__(self, parent, controller, imageframe, plotsframe):
        tk.Frame.__init__(self, parent)
        
        self.controller=controller
        self.imageframe=imageframe
        self.plotsframe=plotsframe
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        
        self.open_file_button=tk.Button(self, text='Open File', 
                                        command = self.controller.open_file, height=2, width=20)
        self.open_file_button.pack(side=tk.TOP, padx=30, pady=10)
        

        
        
        self.cancel_data_processing_button=tk.Button(self, text='Cancel Data Processing', 
                                        command= self.imageframe.cancel_data_processing, height=2, width=20)
        self.cancel_data_processing_button.pack(side=tk.TOP, padx=30, pady=10)
        
        
        self.show_image_button=tk.Button(self, text='Show Data', 
                                        command= lambda: self.controller.show_frame(ACSImageFrame), height=2, width=20)
        self.show_image_button.pack(side=tk.TOP, padx=30, pady=10)
        
        self.show_plots_button=tk.Button(self, text='Show Plots', 
                                        command= lambda: self.controller.show_frame(ACSPlotsFrame), height=2, width=20)
        self.show_plots_button.pack(side=tk.TOP, padx=30, pady=10)

        

        
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



class ACSImageFrame(tk.Frame):
    
    def __init__(self, parent, controller):

        tk.Frame.__init__(self, parent)
        self.controller=controller
        self.timeIn = []
        self.dataArray = []
        self.imagesArray = []
        
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
        
    
    def set_control_frame(self, control_frame):
        self.control_frame=control_frame


    def plot_data(self, image):
        self.a.clear()
        self.im = self.a.imshow(image)
        self.canvas.draw()


    def update_framenum(self, framenum):
        self.current_framenum = framenum
        self.framenum_label.config(text='Frame: '+str(self.current_framenum+1)+'/'+str(self.total_framenum))


    def update_filepath(self, filepath):
        self.filepath = filepath
        self.filepath_label.config(text='Filepath: '+self.filepath)


    def initialize_new_data(self, image, numframes, filepath):
        self.plot_data(image)
        self.total_framenum = numframes
        self.update_framenum(self.current_framenum)
        self.update_filepath(filepath)



    def cancel_data_processing(self):
        self.process_data=False
        self.read_data=False



class ACSPlotsFrame(tk.Frame):
    
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        
        self.controller=controller
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
        
    
    
    def set_control_frame(self, controlframe):
        self.controlframe=controlframe
        return None
            
    def set_slope_struct_functions(self, slope_struct_funcs, r0_vals):
        
        self.slope_struct_functions=slope_struct_funcs
        self.r0=r0_vals
    
    def plot_data(self):
        self.a.clear()
        d_sub=self.controller.frames[ACSDataFrame].pixel_size*self.controller.frames[ACSDataFrame].subaperture_size*self.controller.frames[ACSDataFrame].magnification
        
        rs=np.arange(len(self.slope_struct_functions[self.framenum]))
        r_form=np.linspace(0,len(self.slope_struct_functions[self.framenum]),10000)
        
        D_form=[0.77*6.88/(d_sub**2)*(d_sub/(self.r0[self.r0_framenum]))**(5/3)*(r**2/(1+r**2)+0.438*r**(0.2555)) for r in r_form]
        
        self.x=self.a.plot(rs, self.slope_struct_functions[self.framenum], label=r'$D$')
        
        self.x=self.a.plot(r_form, D_form, label=r'$D_{form}$')
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
        if(self.framenum<len(self.slope_struct_functions)-1):
            self.framenum+=1
            self.r0_framenum+=1
            self.plot_data()


if __name__=='__main__':
    app = ACSDataApp()
    app.mainloop()