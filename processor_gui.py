import tkinter as tk
from tkinter.filedialog import askopenfilename, asksaveasfilename, askdirectory
import os

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
        self.frames[ACSDataFrame].cancel_data_processing()
        time.sleep(0.3)
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
        controls_frame=ACSControlsFrame(container, self, clickerframe, plotsframe)
        
        
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
        
        self.protocol("WM_DELETE_WINDOW", app.safe_quit)
        
        self.geometry('1250x750')

    
    def show_frame(self, cont):
        frame=self.frames[cont]
        frame.tkraise()



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
        
        

        self.numframes_averaged_label=tk.Label(self, text='Number of Frames to Average')
        self.numframes_averaged_label.pack(side=tk.TOP, padx=10, pady=(10,0))
        self.numframes_averaged_entry=tk.Entry(self)
        self.numframes_averaged_entry.insert(0, '60')
        self.numframes_averaged_entry.pack(side=tk.TOP, padx=10, pady=0)

        
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
        
        
        
        self.num_framesprocessed=0
        self.num_framesread=0
        
    
    def set_control_frame(self, controlframe):
        self.controlframe=controlframe
        return None


    def thread_process_data(self):
        threading.Thread(target=self.process_data_method).start()
    

    def process_data_method(self):
        self.process_data=True
        if(self.process_data):
            self.status_label.config(text='Status: Calculating centroids...')
            self.calculate_centroids()
        if(self.process_data):
            self.calculate_cn2()
        # if(self.process_data):
        #     self.status_label.config(text='Status: Plotting data...')
        #     self.controller.frames[ACSPlotsFrame].plot_data()
        #     self.controller.frames[ACSSlopeStructPlotFrame].set_slope_struct_functions(self.slope_struct_functions, self.r0_calculated)
        #     self.controller.frames[ACSSlopeStructPlotFrame].plot_data()
        if(self.process_data):
            self.status_label.config(text='Status: Processing Complete')
        else:
            self.status_label.config(text='Status: Processing Cancelled')
    
    
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