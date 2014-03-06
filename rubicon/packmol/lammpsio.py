# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 22:40:29 2013

@author: asharma6
"""

import re
import numpy as np


class LammpsLog:
    """Parser for LAMMPS log file (parse function). Saves the output properties (log file) in the form of a dictionary (LOG) with the key being
    the LAMMPS output property (see 'thermo_style custom' command in the LAMMPS documentation).
    For example, LOG['temp'] will return the temperature data array in the log file.
    """
    def __init__(self, filename, *properties):
        """
        Args:
            filename:
                Filename of the LAMMPS logfile.
            properties:
                Properties you want to be evaluated. Right now, pymatgen only supports viscosity (argument: viscosity) evaluation.
        """
        
        self.filename=filename
        self.properties=properties
        #print properties
        self.LOG={} # Dictionary LOG has all the output property data as numpy 1D arrays with the property name as the key
        self.header=0
#    @staticmethod
    def list2float(self,seq):
        for x in seq:
          try:
            yield float(x)
          except ValueError:
            yield x
            
    def autocorrelate (a):
        b=np.concatenate((a,np.zeros(len(a))),axis=1)
        c= np.fft.ifft(np.fft.fft(b)*np.conjugate(np.fft.fft(b))).real
        d=c[:len(c)/2]
        d=d/(np.array(range(len(a)))+1)[::-1]
        return d

        
    def parselog(self):
        """
        Parses the log file. rajectory file will be added 
                later with the addition of properties based on trajectory analysis.
        """
        minimization = 0 #To avoid reading the minimization data steps
        
        with open (self.filename,'r') as logfile:
            self.total_lines=len(logfile.readlines())
        logfile.close 
        logfile=open(self.filename,'r')
        for line in logfile:
            if 'viscosity' in self.properties:
                cut=re.search('variable\s+viscosity_cut\s+equal\s+([0-9]+)', line)
                if cut: self.cutoff = float(cut.group(1))
                sep=re.search('variable viscosity_separation equal ([0-9]+)', line)
                if sep: self.separation = float(sep.group(1))
                win=re.search('variable viscosity_window equal ([0-9]+)', line)
                if win: self.window = float(win.group(1))
                t=re.search('variable T equal ([0-9]+)', line)
                if t: self.temp = float(t.group(1))
            
            run=re.search('variable\s+run\s+equal\s+([0-9]+)', line)
            if run: self.runlength = float(run.group(1))
            step=re.search('variable timestep equal ([0-9]+)', line)
            if step: self.timestep=float(step.group(1))
            format = re.search('thermo_style.+', line)
            if format: 
              data_format = format.group().split()[2:] 
              total_properties=len(data_format)
            if ''.join(line.split()).lower()=='run${run}': minimization=1
            if all(isinstance(x,float) for x in list(self.list2float(line.split()))) and minimization==1 and len(line.split())==total_properties: break
            self.header=self.header+1
        if 'viscosity' in self.properties:  
            self.nl0= self.cutoff/self.separation+1 
            self.wline=self.window/self.separation+1 #no. of lines of the window
            self.nline=self.runlength/self.separation+1 
        
        rawdata=np.genfromtxt(fname=self.filename, dtype=float, skip_header=int(self.header+self.nl0), skip_footer=int(self.total_lines-(self.header+self.nline)))
        
        for column, property in enumerate (data_format):
            self.LOG[property]= rawdata[:,column]        
        


    
                    
        
        
         
    
    