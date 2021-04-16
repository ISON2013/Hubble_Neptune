# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 14:46:03 2021

@author: Jamie MacQuillin
"""
#Import all the necessary libraries
import matplotlib.pyplot as plt
import numpy as np
from photutils import centroid_com
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename

from tools import *
from tools.mapping import *

"""

"""
def open_img(path, ext = 'sci'):
    #shield against non-string inputs
    path = str(path)
    ext_dict = {
        'sci': 1, 
        }
    if ext in ext_dict:
        ext = ext_dict[ext]
        
    
    file = get_pkg_data_filename(path)
    fitsdata = fits.open(file)
    imgdata = fits.getdata(file, ext=1)
    fitshdr = fitsdata[0].header
    scihdr = fitsdata[1].header
    return imgdata, fitshdr, scihdr   

"""
Functions which aid with the calibration of the images

- AU_to_km: Takes a distance in AU, converts it to km
                n : distance (km)

- calc_IF: Calculates an I/F factor for a non byte-scaled image, defaults
             based on Hubble WFC3/UVIS instrument
                r : distance between sun and planet (km)
         PHOTFLAM : Inverse sensitivity of the instrument erg/cm^-2/ang/electron
        SolarFlux : Flux density of the sun in that wavelength at Earth
        PlateScale: Plate Scale of the detector (arcsec)

- minnaert: Applies a Minnaert Correction to an image of a disc
"""
def AU_to_km(n):
    Dist = n*149598000
    return Dist

def calc_IF(r, PHOTFLAM, SolarFlux, PlateScale = 0.04, Filter = None):
    # Shielding from arrays of length 1
    dist = float(r)
    PHOTFLAM = float(PHOTFLAM)
    SolFlux = float(SolarFlux)
    PlateScale = float(PlateScale)
    
    if Filter == '467':
        #determining I/F Factor
        omega =((0.0495/3600)*(np.pi/180))**2 #3.76E-14
        fractop = ((PHOTFLAM*10)/omega)
        fracbottom = ((AU_to_km(1.0)/AU_to_km(dist))**2)*(SolFlux/np.pi)
        IF = fractop/fracbottom#PHOTFLAM*(10/omega)*(np.pi/SolFlux)*(dist**2)
    elif Filter == '547':
        #determining I/F Factor
        omega =((0.049/3600)*(np.pi/180))**2 #3.76E-14
        fractop = ((PHOTFLAM*10)/omega)
        fracbottom = ((AU_to_km(1.0)/AU_to_km(dist))**2)*(SolFlux/np.pi)
        IF = fractop/fracbottom#PHOTFLAM*(10/omega)*(np.pi/SolFlux)*(dist**2)
    elif Filter == '619':
        #determining I/F Factor
        omega =((0.0497/3600)*(np.pi/180))**2 #3.76E-14
        fractop = ((PHOTFLAM*10)/omega)
        fracbottom = ((AU_to_km(1.0)/AU_to_km(dist))**2)*(SolFlux/np.pi)
        IF = fractop/fracbottom#PHOTFLAM*(10/omega)*(np.pi/SolFlux)*(dist**2)
    else:
        #determining I/F Factor
        omega =((PlateScale/3600)*(np.pi/180))**2 #3.76E-14
        fractop = ((PHOTFLAM*10)/omega)
        fracbottom = ((AU_to_km(1.0)/AU_to_km(dist))**2)*(SolFlux/np.pi)
        IF = fractop/fracbottom#PHOTFLAM*(10/omega)*(np.pi/SolFlux)*(dist**2)
    return IF, PlateScale

def minnaert(img, cent, k, eph, ex_ang=70):
    """
    Function which applies a Minnaert Function Correction to an image, removing limb darkening
    Inputs are:
    - img: the cropped image of the planetary disc
    - cent: the planicentre and radius of the disc, as given by the get_disc routine from Oliver King
    - k: the minnaert parameter to be applied to the image
    - eph: the ephemeris object corresponding to the image
    - ex_ang: the emission angle from which to exclude the outer edges of the disc, default value 70 degrees
    
    Outputs are:
    - img_corr: the Minnaert Corrected image of the disc
    """
    #Shielding from unexpected forms of inputs for k
    k = float(k)
    thetainc = -float(eph['alpha'])
    thetainc %= 360
    
    #find the theta and phi coordinate images for the disc
    x, y, z = img_to_xyz(img.shape, cent[0], cent[1], cent[2])
    theta, phi = xyz_to_thetaphi(x, y, z)
    theta_i, phi_i = xyz_to_thetaphi(x,y,z,theta=thetainc)
    
    #get the mu parameter image
    mu = np.cos(theta*(np.pi/180))
    mu0 = np.cos(theta_i*(np.pi/180))
    
    #apply the correction to the image
    img_corr = img/((mu**(k-1))*(mu0**k))
    img_corr[theta>ex_ang] = np.nan #trims the extreme edges of the planetary disc off
    return img_corr

"""
Functions which aid with navigation of the images

- circ: Cuts out a circle around the image
            img : the full image of the planet
            rad : the radius of the disc of the planet
            
- get_planicentre: Finds the planicentre of a given disc (starred inputs are required)
            img* : the full image of the planet
            rad* : the radius of the disc of the planet
            This routine can iterate if it does not find the planicentre after a default run through, to find the actual disk
            It does this by cropping the image around the first guess, and then running the guess process iteratively until it converges
            iterations : the number of times the iterative process should be run
            croppix : the starting number of pixels each side of the guess at the planicentre that should be left at iteration 1
            Debug : a boolean variable determining if the images are plotted at each iteration, used when setting iteration and croppix
        returns:
            DiscImg : cropped image of the disc (subarray of the imported img)
            planicentre2 : an object containing the 2 coordinates of the centre
                           and the radius of the disc
"""
def circ(img, rad):
    theta = np.linspace(0, 2*np.pi, 100)
    r = np.sqrt(rad**2)
    y = r*np.cos(theta)
    x = r*np.sin(theta)
    return(img, y, x)

def get_planicentre(img, rad, iterations=1, croppix=400, Debug=False,endcrop=10):
    croppedimg = img
        
    for i in range(iterations):
        imgcent = centroid_com(croppedimg)
        
        if Debug:
            print(imgcent)
            plt.figure(figsize=(5,5))
            im5 = plt.imshow(croppedimg, cmap=col)
            plt.axhline(imgcent[0],color=xcol)
            plt.axvline(imgcent[1],color=xcol)
            plt.title('Partial Image '+str(i))
            plt.colorbar(im5)
            plt.show() 
        
        croppedimg = croppedimg[(int(imgcent[1]-croppix)):(int(imgcent[1]+croppix)),(int(imgcent[0]-croppix)):(int(imgcent[0]+croppix))]
        croppedimg = image.exp_despike(croppedimg, nsigma=1)
        croppix = croppix/3
       
    planicentre1 = centroid_com(croppedimg)
    DiscImg = croppedimg[(int(planicentre1[1])-(int(rad)+endcrop)):(int(planicentre1[1])+(int(rad)+endcrop)),(int(planicentre1[0])-(int(rad)+endcrop)):int(planicentre1[0])+((int(rad)+endcrop))]
    planicentre2 = get_disc(DiscImg, r0=rad)
    return DiscImg, planicentre2
