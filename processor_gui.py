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
        
        self.protocol("WM_DELETE_WINDOW", app.safe_quit)
        
        self.geometry('1250x750')
        
        
        
    
    def show_frame(self, cont):
        frame=self.frames[cont]
        frame.tkraise()


    