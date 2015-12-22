"""
    Partition Data Main GUI
    Used to NEE into GPP and ER
    
    Nighttime ER:  determined from exponential relationship between temperature (nocturnal mean) and Fc (nocturnal sum) on nights without missing observations; binned in soil moisture classes
    Daytime ER:  determined from light response curve between Fsi and Fc; binned in temperature classes
    
    """

import sys

sys.path.append('scripts')

import ast
import copy
import datetime
import logging
import numpy
import time
import Tkinter
import pio
import pls
import putils


class qcgui(Tkinter.Frame):
    def __init__(self, master=None):
        Tkinter.Frame.__init__(self, master)
        self.grid()
        self.createWidgets()

    def createWidgets(self):
        self.process2Label = Tkinter.Label(self,text='L3: Night ER')
        self.process2Label.grid(row=0,column=1,columnspan=1)
        self.process3Label = Tkinter.Label(self,text='L4: Day ER & GPP')
        self.process3Label.grid(row=0,column=2,columnspan=1)
        
        self.fileloadLabel = Tkinter.Label(self,text='Xcel ->')
        self.fileloadLabel.grid(row=1,column=0,columnspan=1)
        self.doxl2nc1Button = Tkinter.Button (self, text="Load Corrected Data", command=self.do_xl2ncCall )
        self.doxl2nc1Button.grid(row=1,column=1,columnspan=1)
        self.doxl2nc2Button = Tkinter.Button (self, text="Load Gapfilled Data", command=self.do_xl2ncCall )
        self.doxl2nc2Button.grid(row=1,column=2,columnspan=1)
        
        self.initiateLabel = Tkinter.Label(self,text='Process Data')
        self.initiateLabel.grid(row=2,column=0,columnspan=1)
        self.doL2Button = Tkinter.Button (self, text="Annual Night ER", command=self.do_l3qc )
        self.doL2Button.grid(row=2,column=1,columnspan=1)
        self.doL3Button = Tkinter.Button (self, text="Annual Day ER & GPP", command=self.do_l4qc )
        self.doL3Button.grid(row=2,column=2,columnspan=1)
        
        self.filesave2Label = Tkinter.Label(self,text='-> Xcel')
        self.filesave2Label.grid(row=3,column=0,columnspan=1)
        self.savexL3Button = Tkinter.Button (self, text="Save L3 Day Data", command=self.do_savexL3 )
        self.savexL3Button.grid(row=3,column=1,columnspan=1)
        self.savexL3Button = Tkinter.Button (self, text="Save L4 Data", command=self.do_savexL4 )
        self.savexL3Button.grid(row=3,column=2,columnspan=1)
        
        self.quitButton = Tkinter.Button (self, text="Quit", command=self.do_quit )
        self.quitButton.grid(row=4,column=0,columnspan=1)
        self.progress = Tkinter.Label(self, text='Waiting for input ...')
        self.progress.grid(row=4,column=1,columnspan=2)

    def do_l3qc(self):
        self.cf = pio.loadcontrolfile('controlfiles')
        if len(self.cf)==0:
            self.do_progress(text='Waiting for input ...')
            return
        if putils.cfkeycheck(self.cf,Base='General',ThisOne='InputLevel'):
            InLevel = self.cf['General']['InputLevel']
        else:
            InLevel = 'L3'
        self.ds2 = pio.nc_read_series(self.cf,InLevel)
        self.do_progress(text='Doing partitioning prep...')
        self.ds3 = pls.l3partition(self.cf,self.ds2)
        self.do_progress(text='Finished partitioning prep')
        self.do_progress(text='Saving L3 Partitioning NetCDF data ...')                     # put up the progress message
        if putils.cfkeycheck(self.cf,Base='General',ThisOne='OutputLevel'):
            OutLevel = self.cf['General']['OutputLevel']
        else:
            OutLevel = 'Partitioning'
        pio.nc_write_series(self.cf,self.ds3,OutLevel)                   # save the L3 data
        self.do_progress(text='Finished L3 Partitioning')              # tell the user we are done
        log.info(' Finished saving L3 Partitioning NetCDF data')
    
    def do_l4qc(self):
        self.cf = pio.loadcontrolfile('controlfiles')
        if len(self.cf)==0:
            self.do_progress(text='Waiting for input ...')
            return
        if putils.cfkeycheck(self.cf,Base='General',ThisOne='InputLevel'):
            InLevel = self.cf['General']['InputLevel']
        else:
            InLevel = 'L4'
        self.ds3 = pio.nc_read_series(self.cf,InLevel)
        self.do_progress(text='Partitioning daytime ER and GPP...')
        self.ds4 = pls.l4partition(self.cf,self.ds3)
        self.do_progress(text='Finished partitioning daytime ER and GPP')
        self.do_progress(text='Saving L4 Partitioning NetCDF data ...')                     # put up the progress message
        if putils.cfkeycheck(self.cf,Base='General',ThisOne='OutputLevel'):
            OutLevel = self.cf['General']['OutputLevel']
        else:
            OutLevel = 'Partitioning'
        pio.nc_write_series(self.cf,self.ds4,OutLevel)                   # save the L3 data
        self.do_progress(text='Finished L4 Partitioning')              # tell the user we are done
        log.info(' Finished saving L4 Partitioning NetCDF data')
    
    def do_progress(self,text):
        self.progress.destroy()
        self.progress = Tkinter.Label(self, text=text)
        self.progress.grid(row=4,column=1,columnspan=3)
        self.update()
    
    def do_quit(self):
        self.do_progress(text='Quitting ...')                         # tell the user what we're doing
        log.info(' Quitting ...')
        self.quit()
    
    def do_savexL3(self):
        self.do_progress(text='Exporting L3 Diurnal NetCDF -> Xcel ...')                     # put up the progress message
        if (putils.cfkeycheck(self.cf,'Output','DefaultXl') and self.cf['Output']['DefaultXl'] == 'False'):
            self.cf = pio.loadcontrolfile('controlfiles')
            if len(self.cf)==0:
                self.do_progress(text='Waiting for input ...')
                return
        if putils.cfkeycheck(self.cf,Base='General',ThisOne='InputLevel') and putils.cfkeycheck(self.cf,Base='General',ThisOne='OutputLevel'):
            InLevel = self.cf['General']['InputLevel']
            OutLevel = self.cf['General']['OutputLevel']
        else:
            InLevel = 'L3'
            OutLevel = 'Partitioning'
        
        pio.autonc2xl(self.cf,InLevel,OutLevel)
        self.do_progress(text='Finished L3 Data Export')              # tell the user we are done
        log.info(' Finished saving L3 data')
    
    def do_savexL4(self):
        self.do_progress(text='Exporting L4 Diurnal NetCDF -> Xcel ...')                     # put up the progress message
        if (putils.cfkeycheck(self.cf,'Output','DefaultXl') and self.cf['Output']['DefaultXl'] == 'False'):
            self.cf = pio.loadcontrolfile('controlfiles')
            if len(self.cf)==0:
                self.do_progress(text='Waiting for input ...')
                return
        if putils.cfkeycheck(self.cf,Base='General',ThisOne='InputLevel') and putils.cfkeycheck(self.cf,Base='General',ThisOne='OutputLevel'):
            InLevel = self.cf['General']['InputLevel']
            OutLevel = self.cf['General']['OutputLevel']
        else:
            InLevel = 'L3'
            OutLevel = 'Partitioning'
        
        pio.autonc2xl(self.cf,InLevel,OutLevel)
        self.do_progress(text='Finished L4 Data Export')              # tell the user we are done
        log.info(' Finished saving L4 data')
    
    def do_xl2ncCall(self):
        self.do_progress(text='Load xl2nc Control File ...')
        self.cf = pio.loadcontrolfile('controlfiles')
        if len(self.cf)==0:
            self.do_progress(text='Waiting for input ...')
            return
        self.do_progress(text='Importing Xcel file -> NetCDF v4 ...')
        if putils.cfkeycheck(self.cf,'General','InLevel') and putils.cfkeycheck(self.cf,'General','OutLevel'):
            InLevel = self.cf['General']['InLevel']
            OutLevel = self.cf['General']['OutLevel']
        else:
            InLevel = 'L3'
            OutLevel = 'L3'
        pio.autoxl2nc(self.cf,InLevel,OutLevel)
        self.do_progress(text='Finished Data Ingest')
        log.info(' Finished Data Ingest')


if __name__ == "__main__":
    log = putils.startlog('partition','logfiles/partition.log')
    qcGUI = qcgui()
    qcGUI.master.title("Carbon Partitioning Main GUI")
    qcGUI.mainloop()
    qcGUI.master.destroy()

    print 'QC: All done'
