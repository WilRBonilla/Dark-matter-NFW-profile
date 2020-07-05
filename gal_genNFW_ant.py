# -*- coding: utf-8 -*-
"""
This code quite simply "generates" a set amount of galaxies seen in the sky
and uses the number of galaxies found to generate random coordinates as well
as their major and minor axes. Finally, it creates a plot of the generated data.

@author: William
"""

# Should be fairly obvious by now, just importing necessary modules.
from pylab import figure, show, rand
from numpy import array, empty, sqrt, arctan, rad2deg, deg2rad, pi, cos, sin, array_str, conjugate, angle
from matplotlib.patches import Ellipse
from numpy.random import normal  
import random 
import numpy as np
from math import atan2


num_galaxies = 1000 # This will be how many galaxies are "found"
# in our patch of sky. Since we do not have actual data, I will use a random, arbitrary value.

print "We have found " + str(num_galaxies) + " galaxies."

galaxy = [] # Empty list to fit however many galaxies we find
sgal = []


counter = 0 # initiate the loop counter as one
index = 0 # initiate the nested loop counter as zero

# The following lines create the empty arrays for each parameter of a galaxy.
# We will save its x & y positions; its major & minor axes; and theta, the angle
# with respect to the x axis.
e1 = empty([num_galaxies,1], dtype ='cfloat')
eg = empty([num_galaxies,1], dtype = 'cfloat')
ef1 = empty([num_galaxies,1], dtype = 'cfloat')
r = empty([num_galaxies, 1], dtype = 'f') # Magnitude of position vector.
g = empty([num_galaxies,1], dtype = 'f') # Magnitude of gamma, or shear effect
R = empty([num_galaxies,1], dtype = 'f') # R is the ratio between minor and major axes
R0 = empty([num_galaxies,1], dtype = 'f')
phi = empty([num_galaxies,1], dtype = 'f')
major = empty([num_galaxies, 1], dtype = 'f')
major0 = empty([num_galaxies, 1], dtype = 'f')
minor = empty([num_galaxies, 1], dtype = 'f')
minor0 = empty([num_galaxies, 1], dtype = 'f')
theta = empty([num_galaxies, 1], dtype = 'f') #angle that galaxy is placed at
theta0 = empty([num_galaxies, 1], dtype = 'f')
eilist = []
#========= Miyoung's Contribution ==============#
Mpc=3.08567758*10**22 #meters
clight=2.99792*10**8/Mpc #m/s -> Mpc/s
z1=0.2 #redshift of lensing halo
z2=1.  #redshift of source
c=5. #concentration
delta_c=(200./3.)*(c**3)/(np.log(1.+c)-c/(1.+c))
Ds= 1.6518*1000. #Gpc at z=1   into Mpc
Dd= 0.6806*1000. #Gpc at z=0.2 into Mpc
Dds= (Ds*(1.+z2)-Dd*(1.+z1))/(1.+z2)
h=0.68 #no unit
H0=1000.*100.*h/Mpc # /s
G=6.67428*10**(-11)/(Mpc**3) # m3/(kg·s^2) --> Mpc^3/(kg s^2)
Sigma_crit= (clight**2)/(4.*np.pi*G)*(Ds/(Dd*Dds)) #kg/Mpc^2
rho_crit = 3.*(H0*H0)/(8.*np.pi*G)*((1.+z1)**3) #Omega_m=1, Omega_lambda=0?
rho_s= rho_crit*delta_c
#M200= (10**14.5) #*(1.989*10**30) #Solarmass into kg
#r200= 1.63*10**(-2)*(M200*h)**(1./3.)/(1.+z1)/h/1000.#kpc -> Mpc
r200= 1.5 #Mpc
M200= (4./3.)*np.pi*r200**3*(200.*rho_crit)
# Assumed Omega(0.2)=1
r_s= r200/c #r_scale

def Sigma_nfw(x):
    if x<1:
        return (1.-(2./np.sqrt(1.-x*x))*np.arctanh(np.sqrt((1.-x)/(1.+x))))/(x*x-1.)
    elif x==1:
        return 1./3.
    elif x>1:
        return (1.-(2./np.sqrt(x*x-1.))*np.arctan(np.sqrt((x-1.)/(x+1.))))/(x*x-1.)
def Sigma_bar_nfw(x):
    if x<1:
        return (2./np.sqrt(1.-x*x)*np.arctanh(np.sqrt((1.-x)/(1.+x)))+np.log(x/2.))/(x*x)
    elif x==1:
        return (1.+np.log(0.5))
    elif x>1:
        return (2./np.sqrt(x*x-1.)*np.arctan(np.sqrt((x-1.)/(1.+x)))+np.log(x/2.))/(x*x)
fac1= (2.*r_s*delta_c*rho_crit)/Sigma_crit
def kappa_nfw(x):
    return fac1*Sigma_nfw(x)
def gamma_nfw(x): #kappa_bar - kappa
    return 2.*fac1*Sigma_bar_nfw(x) - kappa_nfw(x)
def mue(x):
    return 1./np.absolute((1.-kappa_nfw(x))**2-(gamma_nfw(x)**2))
def Gee(x):
    return gamma_nfw(x)/(1 - kappa_nfw(x))

#==========================#

# Start an empty array that holds two values.
xy = empty([1, 2])
position = [xy  for i in range(num_galaxies)] # Initiate a list of galaxies found,
# for now this is a list full of empty arrays.

# Creates a new file or replaces one that already exists. Because we don't have real data,
# This will do for now. Obviously we wont want to replace a file we get from
# an actual telescope. :P
ff = open('telescope.txt', 'w+') 
graphed = 0
# This loop iterates once for each galaxy found

# Mean and the variance of Gaussian distribution
mu, sigma = 0, 0.3
ei1, ei2 = np.random.normal(mu, sigma, num_galaxies), np.random.normal(mu, sigma, num_galaxies)
# e1 is a list of intrinsic ellipticities

for ind in range(0, num_galaxies):
    e1[ind] = np.complex (ei1[ind], ei2[ind])
    #e1[ind] = np.complex (0, 0)


while counter < num_galaxies:
    # r here is the magnitude of the position vector from the origin.
    # Generate a random number with a minimum value of 3 and a maximum value of 10sqrt(2)
    # We want to exclude anything within a 3 parsecs of the source image. 
    #r[counter]= random.uniform(0,5*r_s)

    
    # Find the x and y positions using the polar coordinates from above.
    # xy = array([r[counter]*cos(phi[counter]), r[counter]*sin(phi[counter])])
    xy = array([random.uniform(-2,2), random.uniform(-2,2)])
    position[counter] = [xy] # Simply saves the position for database purposes.
    r[counter] = sqrt(xy[0]**2 + xy[1]**2)  
    y = r[counter]/(r_s) # Dimensionless parameter
    phi[counter] = atan2(xy[1],xy[0])    
    
    # Magnitude of the reduced shear. Gee is the function to find the reduced shear.
    g[counter] = Gee(y)
    
    # This is the complex shear found using |g|exp(2*phi).
    eg[counter] = complex(-g[counter]*cos(2.*phi[counter]), -g[counter]*sin(2.*phi[counter]))

    #Now that we have both shear and source complex ellipticities, to get the final orientation
    #of the galaxy we simply add together the two complex numbers
 
    ef1[counter] = (e1[counter]+eg[counter])/(1 + conjugate(eg[counter])*e1[counter])
    
    # theta is the ellipse angle. This is here only for the purpose of the Ellipse
    # function used ahead and to make sure the shear is applied tangentially, not radially.
    theta[counter] = .5*angle(complex(ef1[counter].real, ef1[counter].imag), deg = 'true')
    theta0[counter] = .5*angle(complex(e1[counter].real, e1[counter].imag), deg = 'true')
        
    #Calculating the magnitude of the complex number modelling the ellipse
    mag = abs(ef1[counter])
    mag0 = abs(e1[counter])
    #This is the expression for how magnitude and axis ratio relate
    R[counter] = (1.-mag)/(1.+mag)
    R0[counter] = (1.-mag0)/(1.+mag0)
    #This generates the actual values for the axes, setting the major to some random
    #number from 1 to .1   
    a = rand(1)*.05
    
    minor[counter] = a
    major[counter] = minor[counter]/R[counter]
    minor0[counter] = a
    major0[counter] = minor0[counter]/R0[counter]
    if mag <= .8: 
        if r[counter] >= .05:
            # Converts each of the elements to a string
            #s = array_str(xy[0]) + "\t"+array_str(xy[1])+"\t"+array_str(ef1[counter][0].real)+"\t"+array_str(ef1[counter][0].imag)+"\n"
            #s =  + "\t"++"\t"++"\t"+ +"\n"
        
            # Fancy way of writing to file, separated by tabs, right justified, and reducing the decimal points.
            s = "{:>10}\t{:>10}\t{:>10}\t{:>10}\n".format('%.5f'%(xy[0]), '%.5f'%(xy[1]),'%.5f'%(ef1[counter][0].real) ,'%.5f'%(ef1[counter][0].imag))        
            ff.write(s) # Write this to the file
            # sgal and galaxies are the source and lensed galaxies, respectively. These
            # lines create lists of the ellipse objects to plot.
            sgal.append(Ellipse(xy, major0[counter], minor0[counter], theta0[counter]))
            galaxy.append(Ellipse(xy, major[counter], minor[counter], theta[counter]))
            graphed += 1
    counter += 1 # No infinite loops. While loop

ff.close() # 
# Plot the locations of galaxies.
fig = figure()
# Adds the subplot of axes to the figure.
ax = fig.add_subplot(111, aspect='equal')
# This is the mass at the center causing the g-lensing.
mass = Ellipse(array([0,0]),r_s,r_s,0)
# For every element, e, in the list galaxy, iterate once.
print "There were " + str(graphed) + " acceptable galaxies"


# Plot the lensed galaxies    
for e in galaxy:
    # Add the current element to axis.
    ax.add_artist(e)
    # Set the face and edgecolors.
    e.set_facecolor('white')
    e.set_edgecolor('blue')
    
# Plot the source galaxies over the lensed ones.
for w in sgal:
    ax.add_artist(w)
    w.set_facecolor('white')
    w.set_edgecolor('red')

# Set maximum and minimum values for the x and y axes and their corresponding
# labels.
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
ax.set_xlabel('x Megaparsec')
ax.set_ylabel('y Megaparsec')

show()

