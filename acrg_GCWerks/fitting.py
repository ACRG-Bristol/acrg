# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 15:06:16 2015

@author: as13988
"""
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import pdb

# Code to generate data
class generate_data:
    def __init__(self,):
        """
        Generate some data 
        """
        x = np.arange(100)
        #print x
        y = -5*np.exp((-0.1*x))-6
        noise = [np.random.random_sample() for i in x]
        noise = np.asarray(noise)
        noise = noise/10+0.95
        
        x1 = [float(i) for i in x]
        x1 = np.asarray(x1)
        y1 = y*noise
        
        print x1.mean()
        print y1.mean()
    
        self.x = x1
        self.y = y1
        
        
# Code to plot data
class plot_data:
    def __init__(self, x, y, fit=0):
        """
        Plot your data
        """
        plt.plot(x, y, 'ro', label="Original Data")
        
        if sum(fit) != 0:
            plt.plot(x, fit, 'b', label="Data Fit")

        plt.legend(loc='upper left')
        plt.show()
        

# Code to fit data
class fit_data:
    def __init__(self, x, y, scale = 0, fit_type=0, sigma = None):

        """
        brutal force to avoid errors
        """    
        x1 = [float(xn) for xn in x] #every element (xn) in x becomes a float
        y1 = [float(yn) for yn in y] #every element (yn) in y becomes a float
        x1 = np.array(x1) #transform your data in a numpy array, 
        y1 = np.array(y1) #so the curve_fit can work
        
        """
        print x1
        print y1
        """
        
        
        # Scale if required
        if scale != 0:
            x_scale = 1/np.mean(x1)
            y_scale = 1/np.mean(y1)
            
            x1 = x_scale*x1
            y1 = y_scale*y1
            if sigma != None:
                sigma = np.array(sigma)
                sigma = y_scale*sigma
        else:
            x_scale = 0
            y_scale = 0
 
        #print x1
        #print y1

       
        """
        make the curve_fit
        """
        #popt, pcov = curve_fit(func, x1, y1) 
                
        if fit_type == 'lin':
            popt,pcov = curve_fit(linfunc, x1, y1, sigma = sigma) 
            fit_label = 'y = ' + str(popt[0]) + 'x + ' + str(popt[1])
            
        else:        
            try:
                popt,pcov = curve_fit(expfunc, x1, y1, sigma = sigma) 
                            
                fit_type = 'exp'
                fit_label =  'y = ' + str(popt[0]) + '*e^-' + str(popt[1])+ 'x +' + str(popt[2])
                
            except RuntimeError:
                print("Error - exp curve_fit failed.")        
                try :
                    popt,pcov = curve_fit(logfunc, x1, y1, sigma = sigma) 
    
                    fit_type = 'log'
                    fit_label = 'y = ' + str(popt[0]) + '*log(-' + str(popt[1])+ 'x) +' + str(popt[2])
                    
                except RuntimeError:
                
                    print("Error - log fit failed. Using linear fit instead.")        
                    popt,pcov = curve_fit(linfunc, x1, y1, sigma = sigma)
    
                    fit_type = 'lin'
                    fit_label = 'y = ' + str(popt[0]) + 'x + ' + str(popt[1])

        if scale !=0 :
            fit = func(x1, popt, fit_type = fit_type)/y_scale
            #print fit
        
        else:
            fit = func(x1, popt, fit_type = fit_type)
            #print fit



        """
        Print the coefficients and plot the funcion.
        """
        print popt
        
        self.fit = fit
        self.coeffs = popt
        self.x_scale = x_scale
        self.y_scale = y_scale
        self.fit_type = fit_type
        self.fit_label = fit_label
        
"""
create a function to fit with your data. a, b and c are the coefficients
that curve_fit will calculate for you. 
In this part you need to guess and/or use mathematical knowledge to find
a function that resembles your data
"""


def expfunc(X, a, b, c):
    return a*np.exp(-b*(X)) + c

def logfunc(X, a, b, c):
    return a*np.log(b*(X)) + c
    
def linfunc(X, a, b):
    return a*X + b
    
def func(X, coeffs, fit_type = 'exp'):
    
    if fit_type =='exp':
        out = expfunc(X, coeffs[0], coeffs[1], coeffs[2])
    
    if fit_type =='log':
        out = logfunc(X, coeffs[0], coeffs[1], coeffs[2])

    if fit_type =='lin':
        out = linfunc(X, coeffs[0], coeffs[1])

    return out