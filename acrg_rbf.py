# -*- coding: utf-8 -*-
"""
Created on Wed May 20 17:05:10 2015

@author: as13988
"""

from scipy.interpolate import Rbf as interp_rbf
import acrg_MOZART as MOZART
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import pdb
import math, random
from copy import copy
from acrg_grid import haversine
from mpl_toolkits.basemap import Basemap 
import datetime

# Convert lat and lons to x, y and z coords
# R is the radius of the earth in km
# Inputs for np.cos, np.sin etc need to be in radians
# Code assumes that the lats/lons are given in degrees
class calc_cartesian():
    def __init__(self, lat, lon, R = 6371):
    
        print 'Assuming that the inputs are in degrees!'
        
        lat = np.radians(lat)
        lon = np.radians(lon)
            
        x = R * np.cos(lat) * np.cos(lon)
        y = R * np.cos(lat) * np.sin(lon)
        z = R * np.sin(lat)
        
        self.x = x
        self.y = y
        self.z = z

# Convert x, y and z coords to lat, lons and radius
# R is in km
class calc_polar():
    def __init__(self, x, y, z):
    
        R = np.sqrt(x**2 + y**2 + z**2)
        lon = np.arccos(z/(R))
        lat = np.arctan(y/x)
        
        self.R = R
        self.lon = lon
        self.lat = lat

# Calculates evenly distributed no. of lats/lons on a sphere
def fibonacci(N, plot = None): 
    inc = np.pi * (3 - np.sqrt(5)) 
    off = 2. / N 
    r2d = 180./np.pi 
    k = np.arange(0,N) 
    y = k*off - 1. + 0.5*off 
    r = np.sqrt(1 - y*y) 
    phi = k * inc 
    x = np.cos(phi)*r 
    z = np.sin(phi)*r 
    theta = np.arctan2(np.sqrt(x**2+y**2),z) 
    phi = np.arctan2(y,x) 
    lats = 90.-r2d*theta 
    lons = r2d*phi 
    
    if plot != None:
        map = Basemap(projection ='ortho',lat_0=0,lon_0=-90) 
        map.drawcoastlines() 
        map.fillcontinents(color='coral') 
        map.drawmapboundary(fill_color='aqua') 
        x,y = map(lons, lats) 
        map.scatter(x,y,10,marker='o',zorder=10) 
        plt.show()
    
    return lats, lons 


# Calculate a set of n points (n=samples) evenly spread around a globe in cartesian coords
class fibonacci_sphere():
    def __init__(self, samples=1,randomize=True, out_filename = None, title = None, plot = None):
        rnd = 1.
        if randomize:
            rnd = random.random() * samples
    
        points = []
        offset = 2./samples
        increment = math.pi * (3. - math.sqrt(5.));
    
        for i in range(samples):
            y = ((i * offset) - 1) + (offset / 2);
            r = math.sqrt(1 - pow(y,2))
           
            
            phi = ((i + rnd) % samples) * increment
    
            x = math.cos(phi) * r
            z = math.sin(phi) * r
    
            points.append([x,y,z])
    
        if plot is not None:
            # plot the output
            fig_size = plt.rcParams["figure.figsize"]
            fig_size[0] = 6
            fig_size[1] = 6
            plt.rcParams["figure.figsize"] = fig_size
            
            fig = plt.figure()
            ax = Axes3D(fig)
    
            
    
            x_s=[];y_s=[]; z_s=[]
        
            for pt in points:
                x_s.append( pt[0]); y_s.append( pt[1]); z_s.append( pt[2])
        
            ax.scatter3D( np.array( x_s), np.array( y_s), np.array( z_s) )                
            ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
            
            if title is not None:
                ax.set_title(title, fontsize=16)
    
            plt.show()
            
            if out_filename is not None:
                fig.savefig('/home/as13988/Plots/' + out_filename)

        # Consolidate the output
        self.x = np.array(points)[:,0]
        self.y = np.array(points)[:,1]
        self.z = np.array(points)[:,2]
        self.points = np.array(points)

class map_E_space():
        def __init__(self, e_0, xi, ci, wi, y, fn_type ='gaussian'):
               
            J = np.zeros(100)   
            E = np.random.randint(0, 50, 100)/10.0
            E[0] = e_0
            
            for i in np.arange(100):
               
                out_i = rbf_fn_E(E[i], xi, ci, wi)
                
                J_i = (y - out_i)**2
                
                J[i] = np.mean(J_i)
           
                print 'E = ' + str(E[i]) + ', J = ' + str(J[i])
                
            self.J = J
            self.E = E

# xi = the values
# ci = the centre
# wi = the weights
# epsilon = the width of the RBF's
# This is just a function wrapper for calc_multi_rbf with the epsilon listed as the main arguement
def rbf_fn_E(epsilon, xi, ci, wi, fn_type ='gaussian'):
    
    out = calc_multi_rbf(xi, ci, wi = wi, epsilon = epsilon, fn_type = 'gaussian')
            
    return  out.sum

# xi = the values
# ci = the centre
# wi = the weights
# epsilon = the width of the RBF's
# This is just a function wrapper for calculating the cost function J with the epsilon listed as the main arguement
def J_fn_E(epsilon_wi, xi, ci, y, fn_type ='gaussian'):

    print 'Called J_fn_E'   

    epsilon = epsilon_wi[0:len(epsilon_wi)/2]
    wi = epsilon_wi[len(epsilon_wi)/2:len(epsilon_wi)]
    
    print 'wi'
    print wi
    print 'epsilon'
    print epsilon
    
    rbf = calc_multi_rbf(xi, ci, wi = wi, epsilon = epsilon, fn_type = 'gaussian')
    
    J = np.mean((y - rbf.sum)**2)
    
    print 'J = ' + str(J)
    
    #pdb.set_trace()    
    
    return  J

# minimise the cost (y  - yi)**2 based on E using optimize.minimize
# E_0 = the inital guess at Epsilon
# xi = the inputs
# ci = the locations of the centres
# wi = the weighting of each RBF
# y = the output that we're aiming for

class rbf_E_optimise():
    def __init__ (self, E_0, xi, ci, wi, y):
        from scipy import optimize
        
        # Initalise the epsilon value to the inital guess
        out = optimize.minimize(J_fn_E, E_0, args = (xi, ci, wi, y))
        
            
        self.out = out

class rbf_wiE_optimise():
    def __init__ (self, E_0, wi, xi, ci, y):
        from scipy import optimize
        
        # Initalise the epsilon value to the inital guess
        out = optimize.minimize(J_fn_E, [E_0, wi], args = (xi, ci, y), method='Nelder-Mead')
        
            
        self.out = out

# m denotes the number of examples here, not the number of features
# E_0 = the inital guess at Epsilon
# xi = the inputs
# ci = the locations of the centres
# wi = the weighting of each RBF
# y = the output that we're aiming for
# cost_bound = the size of the mean difference at which we'll call it done
# diff_bound = the point a which the difference between iterations of epsilon becomes negligible
class rbf_gradientdescent():
    def __init__ (self, E_0, xi, ci, wi, y, lat, lon, eta, Iterations, cost_bound = 0.1, diff_bound = 0.1, plot = None, quiet =None):
    
        # reform the output from lat*lon to an array of 1D
        y = y.reshape(y.size)    
    
        # Step 1 - Initalise the epsilon value to the inital guess
        E_next = np.zeros(len(wi))
        E_next[:] = E_0
        
        print 'This code assumes a gaussian RBF'
        print "If you want to use another funtion type then you'll need to rewrite the code."
        
        cost = cost_bound
        diff = diff_bound
        
        E_array = np.zeros((Iterations, len(wi)))
        gradients = np.zeros((Iterations, len(wi)))
        cost_array = np.zeros(Iterations)
         
            
        j = 0        
        
        while cost >= cost_bound and np.max(diff) >= diff_bound and j <= Iterations-1:
            # Step 2 - Calculate the output value using the given E_i, xi, ci and wi
            # F(E) = w * exp (-(Er)^2)
            E_j = copy(E_next)
            E_array[j,:] = E_j 
            
            rbf_j = calc_multi_rbf(xi, ci, wi = wi, epsilon=E_j)
                
            yj = rbf_j.sum
            
            # plot the input and the rbf out for for comparison
            if plot != None:
                lon_gr, lat_gr = np.meshgrid(lon, lat) 
                
                plot_surf(lon_gr, lat_gr, yj.reshape((96,144)), title = 'Fit - j = ' + str(j))
            
        
            
            # Step 3 - Calculate the "cost" function and the average cost
            J = np.mean((y - yj)**2)
             
            cost = np.sum(np.sqrt(J))/J.size      #  The average cost over the whole field          
            cost_array[j] = cost
            
            print("Iteration %d | Cost: %f" % (j, cost))
            print("Iteration %d | J: %f" % (j, np.sum(J)))
                        
            
            # Step 4 - Calculate the gradient
            # J = (y - yi)**2
            # J = (y - w * exp (-*E*r)**2) )**2
            # dJ/dE = 4*E*(y -yi).[w.((r^^2)exp (-(E*r)**2))]
            #pdb.set_trace()
            #gradient = 4*E_j*np.dot((y-yj), np.dot(wi, np.multiply((rbf_j.r**2),rbf_j.phi)))
            # Calculate the gradient for each E_j
            dist = y - yj            
            
                        
            
            for i in np.arange(len(E_j)):
                
                rji = rbf_j.r[i,:]
                rji_sqrd = rji**2
                phi = rbf_j.phi[i,:]             
                
                
                gradient_j = 4*E_j[i]*np.dot(np.transpose(dist),wi[i]*np.multiply(rji_sqrd,phi))
                gradients[j,i] = gradient_j
                
                # Calculate the next value of E
                E_next[i] = E_j[i] - eta * gradient_j
            
            pdb.set_trace()        
            
            # Calculate the difference between E_i values            
            diff = np.abs(E_next - E_j)
            
            j = j + 1
            
            print 'E max diff = ' + str(np.max(diff))
            
            print 'E next = ' +str(E_next)
                
                
        
        print  'Cost = ' + str(cost >= cost_bound)
        print 'E max diff = ' + str(np.max(diff) >= diff_bound)
        print 'E max diff = '+ str(np.max(diff))
        print 'Iterations = ' + str(j <= Iterations-1)
            
        self.E = E_j
        self.E_array = E_array
        self.cost_array = cost_array
        self.gradients = gradients
        self.fit = calc_multi_rbf(xi, ci, wi = wi, epsilon=E_j).sum

# Forth case
# Uses the code from the third case to:
#    - fix the number of centres then evenly distribute around the globe
#    - calculate a first estimate using an inital fixed RBF width (Epsilon) for each centre
#    - calculate an inital first guess for the weights
# lats/lons converted to radial coords
# It then alternates between:
#   1) using gradient descent to calculate the best multiple epsilon values for the given wieghts
#   2) recalculates the weights bases on the new epsilon values
class rbf_4():
    def __init__(self, lats, lons, nocentres, output, eta = 0.001, Iterations = 100, cost_bound = 0.1, \
                        diff_bound = 1, epsilon = 0.1, plot = None, no_norm = None):
    
        # Set up the initial calculation of the weights  
        
        # Convert lats and lons to grid
        lon_gr, lat_gr = np.meshgrid(lons, lats)  
        
        # Convert inputs to a vector of pairs
        xi = np.array((lat_gr.flatten(), lon_gr.flatten()))            
        
        # Evenly distribute the centres around the globe in lats/lons
        ci = np.array(fibonacci(nocentres, plot=1))

        # Convert outputs to vector
        y = output.flatten()
        
        yi = y - min(y)         
        
        # Want to start with E = average distance between centres
        # Calculate the distance between the each centre and the rest
        for i in np.arange(36):
            temp = np.zeros(len(ci[0,:]))
            
            # syntax (origin, lat, lon)
            ci_dist = haversine.fn_multipledistances(ci[:,i], ci[0,:], ci[1,:], temp)
            
            if i == 0:
                dist = ci_dist
            else:
                dist = np.concatenate((dist,ci_dist))
        
            
        
        epsilon_i = np.zeros(nocentres)
        epsilon_i[:] = (1/np.mean(dist))**2
        
        
        i = 0
        
        J_w = 0
        J_E = 2
             
        
        while abs(J_w - J_E) >= 1 and i <= Iterations-1:

            # Now need to solve
            # yi = wi phi(xi, ci)
            # where wi are the weights of each phi RBF
            
            # Calculate phi(xi,ci) 
            rbf_i = calc_multi_rbf(xi, ci, epsilon=epsilon_i)
            
            #plot_surf(lon_gr, lat_gr, rbf_i.sum.reshape(len(lats), len(lons)))
            
            phii = rbf_i.phi
            
            #plt.plot(np.arange(13824), phii[0,:])            
            #plt.show()
            
            # we want to solve
            # wi = phii+ yi
            # where phii+ is the pusedoinverse of phii       
            phii_pinv = np.linalg.pinv(phii)
            
            # Calculate the weights from the pseudoinverse
            wi = np.dot(yi,phii_pinv)
            
            # Calculate the fit       
            fit_w = np.dot(wi,phii) + min(y)
            
           
        
            # plot the input and the rbf out for for comparison
            if plot != None:
                if isinstance(epsilon, np.ndarray):
                    epsilon_tag = copy(epsilon[0])
                else:
                    epsilon_tag =copy(epsilon)
                
                if i == 0:                    
                    plot_surf(lon_gr, lat_gr, output, title='Input', out_filename='Fits/RBF_4_input.png')
                
                plot_surf(lon_gr, lat_gr, fit_w.reshape(len(lats), len(lons)), title = 'Fit - Epsilon = ' + str(epsilon_tag), out_filename='Fits/RBF_4_fit_E_' + str(epsilon_tag) + '.png')
        
            #pdb.set_trace()
             
            J_w = np.mean((y - fit_w)**2)
            
            B = copy(epsilon_i)            
            wi_in = copy(wi)            
            
            A = rbf_wiE_optimise(B, wi_in, xi, ci, yi)
            
            pdb.set_trace()
            
            """
            pdb.set_trace()
            
            # Do the gradient descent fit for the Ei's    
            GD = rbf_gradientdescent(epsilon_i, xi, ci, wi, yi, lats, lons, eta, Iterations, \
                cost_bound = cost_bound, diff_bound = diff_bound, plot = plot, quiet = None)
    
            pdb.set_trace()    
    
            fit_GD = GD.fit + min(y)    
    
            J_E = np.mean((y - fit_GD)**2)
        
            epsilon_i = GD.E
            
            
        
            # plot the input and the rbf out for for comparison
            if plot != None:
                if isinstance(epsilon, np.ndarray):
                    epsilon_tag = copy(epsilon[0])
                else:
                    epsilon_tag =copy(epsilon)
                
                if i == 0:                    
                    plot_surf(lon_gr, lat_gr, output, title='Input', out_filename='Fits/RBF_4_input.png')
                
                plot_surf(lon_gr, lat_gr, fit.reshape(len(lats), len(lons)), title = 'Fit - Epsilon = ' + str(epsilon_tag), out_filename='Fits/RBF_4_fit_E_' + str(epsilon_tag) + '.png')
                    
            
            
            i = i + 1
            


        self.GD = GD
        """
        self.A

# Third case
# Fix the number of centres then evenly distribute around the globe
# Fixed RBF width (Epsilon)
# Only fitting weights
# lats/lons converted to radial coords
class rbf_3():
    def __init__(self, lats, lons, nocentres, output, epsilon = 0.1, plot = None, no_norm = None):

        # Convert lats and lons to cartesian co-ordiantes
        lon_gr, lat_gr = np.meshgrid(lons, lats)  
        
        # Use a radius of 1 to match the generated centres
        coords = calc_cartesian(lat_gr, lon_gr, R=1)
        
        # Evenly distribute the centres around the globe
        centre_coords = fibonacci_sphere(samples = nocentres, randomize=False, \
            out_filename = 'Fits/Centres_N_eq_' + str(nocentres) + '.png', \
            title = str(nocentres) + ' centres evenly spaced', plot = 1)  
        

        # Convert inputs to a vector of pairs
        xi = prep(coords.x, input2 = coords.y, input3 = coords.z).output
        ci = np.transpose(centre_coords.points)
        
        # Convert outputs to vector
        yi = prep(output).output
        
        # Renormalise it so that the lowest conc = 0
        min_yi = np.min(yi)
        
        if no_norm is None:
            print 'Deducting min yi value during calculations'
            yi = yi - min_yi
        else :
            yi = yi
            min_yi = 0
        
        # Now need to solve
        # yi = wi phi(xi, ci)
        # where wi are the weights of each phi RBF
        
        # Calculate phi(xi,ci) 
        rbf = calc_multi_rbf(xi, ci, epsilon=epsilon)
        phii = rbf.phi
        
        # we want to solve
        # wi = phii+ yi
        # where phii+ is the pusedoinverse of phii       
        phii_pinv = np.linalg.pinv(phii)
        
        # Calculate the weights from the pseudoinverse
        wi = np.dot(yi,phii_pinv)
        
        # Calculate the fit and add back the minimum concentration                
        fit = np.dot(wi,phii) + min_yi
    
        # plot the input and the rbf out for for comparison
        if plot != None:
            #pdb.set_trace()
            if isinstance(epsilon, np.ndarray):
                epsilon = epsilon[0]
            plot_surf(lon_gr, lat_gr, output, title='Input', out_filename='Fits/RBF_3_input.png')
            plot_surf(lon_gr, lat_gr, fit.reshape(len(lats), len(lons)), title = 'Fit - Epsilon = ' + str(epsilon), out_filename='Fits/RBF_3_fit_E_' + str(epsilon) + '.png')
    
        self.weights = wi
        self.phi = phii
        self.phi_inv = phii_pinv
        self.fit = fit.reshape(len(lats),len(lons))
        self.ci = ci
        self.xi = xi
        self.yi = yi + min_yi
        self.ri = rbf.r
        self.y_fit = fit 

# Second case
# Fixed centres, fixed width
# Only fitting weights
# lats/lons converted to radial coords
class rbf_2():
    def __init__(self, lats, lons, centre_lats, centre_lons, output, epsilon = 0.1, plot = None):

        # Convert lats and lons to cartesian co-ordiantes
        lon_gr, lat_gr = np.meshgrid(lons, lats)  
        coords = calc_cartesian(lat_gr, lon_gr)
        centre_coords = calc_cartesian(centre_lats, centre_lons)


        # Convert inputs to a vector of pairs
        xi = prep(coords.x, input2 = coords.y, input3 = coords.z).output
        ci = prep(centre_coords.x, input2 = centre_coords.y, input3 = centre_coords.z).output
        
        # Convert outputs to vector
        yi = prep(output).output
        
        # Renormalise it so that the lowest conc = 0
        min_yi = np.min(yi)
        
        yi = yi - min_yi
        
        # Now need to solve
        # yi = wi phi(xi, ci)
        # where wi are the weights of each phi RBF
        
        # Calculate phi(xi,ci) 
        phii = (calc_multi_rbf(xi, ci, epsilon=epsilon)).phi
        
        # we want to solve
        # wi = phii+ yi
        # where phii+ is the pusedoinverse of phii       
        phii_pinv = np.linalg.pinv(phii)
        
        # Calculate the weights from the pseudoinverse
        wi = np.dot(yi,phii_pinv)
        
        # Calculate the fit and add back the minimum concentration
        fit = np.dot(wi,phii) + min_yi
    
        # plot the input and the rbf out for for comparison
        if plot != None:
            plot_surf(lon_gr, lat_gr, output, title='Input', out_filename='Fits/RBF_3_input.png')
            plot_surf(lon_gr, lat_gr, fit.reshape(len(lats), len(lons)), title = 'Fit - Epsilon = ' + str(epsilon), out_filename='Fits/RBF_2_fit_E_' + str(epsilon) + '.png')
    
        self.weights = wi
        self.phi = phii
        self.phi_inv = phii_pinv
        self.fit = fit.reshape(len(lats),len(lons))
        self.ci = ci
        self.xi = xi
        


# Simplest case
# Fixed centres, fixed width
# Only fitting weights
# Using lats/lons not converted to radial coords
class rbf_1():
    def __init__(self, lats, lons, centre_lats, centre_lons, output, epsilon = 0.1, plot = None):

        # Convert inputs to a vector of pairs
        xi = prep(lons, input2 = lats).output
        ci = prep(centre_lons, input2 = centre_lats).output
        
        # Convert outputs to vector
        yi = prep(output).output
        
        # Renormalise it so that the lowest conc = 0
        min_yi = np.min(yi)
        
        yi = yi - min_yi
        
        # Now need to solve
        # yi = wi phi(xi, ci)
        # where wi are the weights of each phi RBF
        
        # Calculate phi(xi,ci) 
        phii = (calc_multi_rbf(xi, ci, epsilon=epsilon)).phi
        
        # we want to solve
        # wi = phii+ yi
        # where phii+ is the pusedoinverse of phii       
        phii_pinv = np.linalg.pinv(phii)
        
        wi = np.dot(yi,phii_pinv)
        
        fit = np.dot(wi,phii) + min_yi
    
        # plot the input and the rbf out for for comparison
        if plot != None:
            
            lon_gr, lat_gr = np.meshgrid(lons, lats)  
            plot_surf(lon_gr, lat_gr, output, title='Input', out_filename='Fits/RBF_1_input.png')
            plot_surf(lon_gr, lat_gr, fit.reshape(len(lats), len(lons)), title = 'Fit - Epsilon = ' + str(epsilon), out_filename='Fits/RBF_1_fit_E_' + str(epsilon) + '.png')
    
        self.weights = wi
        self.phi = phii
        self.phi_inv = phii_pinv
        self.fit = fit.reshape(len(lats),len(lons))



# Wrapper to get the inputs or outputs into the appropriate shape for rbf
# Takes either 2 inputs which are currently either:
#       - two vectors, one of length = m and one of length = n
#       - two 2D matricies, one of lats (m by n) and lons (m by n)
# and converts them to a grid to match the output array shape  and then to m * n sets of pairs 
# Or a single input of m x n matrix to a vector of length m * n 
# NB: THIS CAN ONLY COPE WITH INPUTS OF 2D 
# I"LL NEED TO ADAPT THIS FOR HIGHER ORDER PROBLEMS
class prep():
    def __init__(self, input1, input2 = None, input3 = None):

        # only a single input is given
        if input2 == None and input3 == None:    
            output = input1.flatten()

        # two inputs are given
        if input2 is not None and input3 is None:
            
            # check if the inputs are of the same dimension
            if len(np.shape(input1)) == len(np.shape(input2)):
                
                # If they're vectors make matricies
                if len(np.shape(input1)) == 1 and len(np.shape(input2)) == 1:
                    input1_matrix, input2_matrix = np.meshgrid(input1, input2)
                
                # They're already matricies
                else:
                    input1_matrix = input1
                    input2_matrix = input2
                
                # combine the matricies to make a vector of pairs
                output = np.array((input1_matrix.flatten(), input2_matrix.flatten()))
                
            else:
                print "If you're providing two inputs then they need to be either:"
                print " - two vectors, one of length = m and one of length n"
                print " - two 2D matricies, one of lats (m by n) and lons (m by n)"
                
                output = None
                
        # three inputs are given
        if input2 is not None and input3 is not None:
            
            # check if the inputs are of the same dimension
            if len(np.shape(input1)) == len(np.shape(input2))  and len(np.shape(input2)) == len(np.shape(input3)):
                
                # If they're vectors make matricies
                if len(np.shape(input1)) == 1 and len(np.shape(input2)) == 1 and len(np.shape(input3)) == 1:
                    input1_matrix, input2_matrix = np.meshgrid(input1, input2)
                    input1_matrix, input3_matrix = np.meshgrid(input1, input3)
                
                # They're already matricies
                else:
                    input1_matrix = input1
                    input2_matrix = input2
                    input3_matrix = input3
                
                # combine the matricies to make a vector of pairs
                output = np.array((input1_matrix.flatten(), input2_matrix.flatten(), input3_matrix.flatten()))
                
            else:
                print "If you're providing three inputs then they need to be either:"
                print " - three vectors of equal length"
                print " - three matricies of the same shape)"
                
                output = None
       
        
        self.output = output
    

# Vector norm
def euclidean_norm(x1, x2):
    return np.sqrt( ((x1 - x2)**2).sum(axis=0) )

# Norm of the difference between each vector in the matrix (x1) of vectors and a vector (x2)
def calc_r(x1, x2):
    
    shape = np.shape(x1)
    
    r = np.zeros(shape[1])
    
    temp = copy(x1)
    temp[0,:] = x2[0]    
    temp[1,:] = x2[1]    
    temp[2,:] = x2[2]    
    
    r = np.sqrt(np.sum(np.abs(x1-temp)**2,0)) 
    
    # Old loop version is 300 times slower than matrix version above
    #for i in np.arange(shape[1]):
           # r[i] = np.linalg.norm(x1[:,i] - x2)
    
    return r



# x = the values
# ci = the centres
class calc_multi_rbf():
    def __init__(self, xi, ci, wi = None, epsilon=0.1, fn_type ='gaussian'):
        
        # x and xi need to be of the same first dimension
        # e.g. if the centres are vectors of dim 3 than the    
        # values for evaluation need to be vectors of dim 3 either a matrix of these i.e. dim = 3 x n x m
        #  or a vector i.e. dim =  3 x n
          
         
        if np.shape(xi)[0] != np.shape(ci)[0]:
            
            print "x and ci need to be of the same first dimension"
            print "e.g. if the centres are vectors of dim 3 than the values for evaluation need to be vectors " \
                + "of dim 3 either a matrix of these i.e. dim = 3 x n x m or a vector i.e. dim =  3 x n"
            
            outarray = 0
            
        else:        
                    
            nocentres = np.shape(ci)[1]        
    
            #pdb.set_trace()
            if wi == None:
                wi = np.zeros(nocentres)
                wi[:] = 1
                
            # Make an array to store the output values at each point for each centre
            # dim = [no centres, dim x]
            outarray = np.zeros(np.append(nocentres, np.shape(xi)[1:]))       
            
            # Make an array to store the output norm value at each point for each centre
            # dim = [no centres, dim x]
            out_r = np.zeros(np.append(nocentres, np.shape(xi)[1:]))       
            
                                 
            # Loop through each centre
            for i in np.arange(nocentres):
                                
                # Allow it to accept multiple values of epsilon
                if isinstance(epsilon, (frozenset, list, set, tuple, np.ndarray)):
                    if len(epsilon) == 1:
                        epsilon_i = epsilon[0]
                    else:
                        epsilon_i = epsilon[i]
                else:
                    epsilon_i = epsilon
                
                # calculate the rbf for the given centre 
                rbf_i = calc_rbf(xi, ci[:,i], fn_type = fn_type, epsilon = epsilon_i)
                
                outarray[i, :] = rbf_i.out
                
                out_r[i,:] = rbf_i.r
            
            
            # multiply by the weights
            out = np.dot(wi, outarray)

            
        self.sum = out
        self.phi = outarray
        self.r = out_r
        self.wi = wi
        

# x = the values
# xi = the centre
class calc_rbf():
    def __init__(self, xi, ci, epsilon=1.0, fn_type ='gaussian'):
    
        # x and xi need to be of the same first dimension
        # e.g. if the centres are vectors of dim 3 than the    
        # values for evaluation need to be vectors of dim 3 either a matrix of these i.e. dim = 3 x n x m
        #  or a vector i.e. dim =  3 x n
        if np.shape(xi)[0] != np.shape(ci)[0]:
            
            print "x and xi need to be of the same first dimension"
            print "e.g. if the centre is a vector of dim 3 than the values for evaluation need to be vectors " \
                + "of dim 3 either a matrix of these i.e. dim = 3 x n x m or a vector i.e. dim =  3 x n"
        
        else:     
                    
            # r is the radius or the norm
            #r = calc_r(x,xi)
            temp = np.zeros(np.shape(xi[1]))      

                   
            
            if len(ci) == 2 :
                r = (haversine.fn_multipledistances(ci, xi[0,:], xi[1,:], temp))**2
            else:
                r = calc_r(xi, ci)
                
            dictionary = {
                        'multiquadric': fn_multiquadric(r,epsilon), \
                        'inverse': fn_inverse(r,epsilon), \
                        'gaussian': fn_gaussian(r,epsilon), \
                        'linear': fn_linear(r,epsilon), \
                        'quintic': fn_quintic(r,epsilon), \
                        'thin_plate': fn_thin_plate(r,epsilon)}           
            
            """
            
            if fn_type == 'multiquadric': 
                 out = np.sqrt((r/epsilon)**2 + 1)
            
            out = 0           
            if fn_type ==  'inverse':
                out = 1.0/np.sqrt((r/epsilon)**2 + 1)
            
            if fn_type == 'gaussian': 
                out = np.exp(-(r*epsilon)**2)
                #out = np.exp((-1)*epsilon*(r**2))
                
            if fn_type == 'linear': 
                out = r
                
            if fn_type == 'cubic': 
                out = r**3
                
            if fn_type == 'quintic': 
                out = r**5
                
            if fn_type == 'thin_plate': 
                out = r**2 * np.log(r)
             """
             
        self.out =  dictionary[fn_type]
        self.r = r

# The RBF function types
def fn_multiquadric(r,epsilon):
    
    out = np.sqrt((r/epsilon)**2 + 1)
    
    return out
    
def fn_inverse(r, epsilon):
    
    out = 1.0/np.sqrt((r/epsilon)**2 + 1)
    
    return out
    
def fn_gaussian(r, epsilon):
    
    out = np.exp(-(r*epsilon)**2)
    
    return out
    
def fn_linear(r, epsilon): 
    
    out = r
    
    return out
    
def fn_cubic(r, epsilon): 
                
    out = r**3
    
    return out
                
def fn_quintic(r, epsilon): 
    
    out = r**5
    
    return out
                
def fn_thin_plate(r, epsilon): 
    
    out = r**2 * np.log(r)
    
    return out



# Do a quick surface plot
def plot_surf(x, y, z, q = None, scale =1, out_filename = None, title = None):
    
    plt.rcParams["figure.figsize"] = [6,4]   
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
     
    
    if q is None:
        surf = ax.plot_surface(x, y, z*scale, rstride=1, cstride=1, cmap=cm.gist_rainbow, \
            linewidth=0, antialiased=False)
        fig.colorbar(surf, shrink=0.5, aspect=5)
         
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

        ax.set_ylabel('Lat')
        #ax.set_ylabel()  
        ax.set_xlabel('Lon')
        

    else:
        q_norm = q/q.max()
        surf = ax.plot_surface(x, y, z*scale, rstride=1, cstride=1, facecolors=cm.gist_rainbow(q_norm), \
            linewidth=0, antialiased=False)
        m = cm.ScalarMappable(cmap=cm.gist_rainbow)
        m.set_array(q_norm)
        plt.colorbar(m)

        
    if title is not None:
        ax.set_title(title, fontsize=16)

    

    plt.show()
    
    if out_filename is not None:
        fig.savefig('/home/as13988/Plots/' + out_filename)
    
    plt.close()