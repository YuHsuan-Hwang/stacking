#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  7 20:50:13 2019

@author: yuhsuan
"""

# =============================================================================
# packages
# =============================================================================
from __future__ import division
import time
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
import random
# =============================================================================
# functions
# =============================================================================
def ReadCatalog(index,catalog):
    data[index] = Table.read(catalog, hdu=1)
    x[index] = data[index][color2]-data[index][color3]
    y[index] = data[index][color1]-data[index][color2]
    return

def Mask_M(index):
    down = -30
    up = 0
    mask_color1 = (data[index][color1]>down) & (data[index][color1]<up)
    mask_color2 = (data[index][color2]>down) & (data[index][color2]<up)
    mask_color3 = (data[index][color3]>down) & (data[index][color3]<up)
    mask_M = mask_color1 & mask_color2 & mask_color3
    return mask_M

def Mask_error(mode,inputup,index):
    up = inputup
    if mode==1:
        mask_V_errorV = (data[index][V_error]<up) & (data[index][V_error]>0)
        #Suprime-Cam:i+band
        mask_ip_error = (data[index][ip_error]<up) & (data[index][ip_error]>0)
        mask_J_error = (data[index][J_error]<up) & (data[index][J_error]>0)
        mask_Ks_error = (data[index][Ks_error]<up) & (data[index][Ks_error]>0)
        return mask_V_errorV | mask_ip_error | mask_J_error | mask_Ks_error
    if mode==2:
        return (data[index][Ks_error]<up) & (data[index][Ks_error]>0)
    if mode==3:
        return (data[index][Ks_mag]<24)
    

def Mask_photoz(index):
    return (data[index][photoz]>0) & (data[index][photoz]<8)

def Mask_class_star(index):
    return (data[index][class_star]==0)

def Mask_classQG(index):
    return (data[index][class_SFG]==0)
def Mask_classSFG(index):
    return (data[index][class_SFG]==1)
def Mask_myclassQG(index):
    return (y[index]>3.1) & (y[index]>3.0*x[index]+1.0)
def Mask_myclassSFG(index):
    return (y[index]<=3.1) | (y[index]<=3.0*x[index]+1.0)

def MaskRegion(inputra,inputdec,inputradius):
    c0 = SkyCoord(ra=data[0][mask[0]][ra[0]], dec=data[0][mask[0]][dec[0]])
    center = SkyCoord(inputra, inputdec, frame='fk5')
    radius = inputradius
    sep = center.separation(c0)
    return sep<=radius

def ReadImage(filename):
    hdu_list = fits.open(filename)
    #hdu_list.info()
    image_data = hdu_list[0].data[0]
    wcs = WCS(hdu_list[0].header).celestial
    return image_data,wcs

def Stacking(data_tmp,inputvmin,inputvmax,inputrms):
    
    final_image = np.zeros((30,30))
    not_use = 0
    print 'len(data_tmp.filled()) ',len(data_tmp.filled())
   
    picnum = 1
    
    rnds = random.sample(range(0, len(data_tmp.filled())), 9) 
    print rnds
    
    for i in range( len(data_tmp.filled()) ):
        ra1 = data_tmp[ra[0]][i]
        dec1 = data_tmp[dec[0]][i]
        x1, y1 = wcs.all_world2pix([[ra1, dec1]], 0)[0]
        image = image_data[int(y1)-15:int(y1)+15:1,int(x1)-15:int(x1)+15:1]
    
        if (image.shape==(30,30)) & (np.isnan(image).any()==False):
            final_image += image
        else:
            not_use += 1
        if i in rnds:
            plt.subplot(3, 3, picnum)
            picnum = picnum+1
            if inputvmax==9999:
                plt.imshow(image, cmap='gray', vmin=inputvmin, vmax=image.max())
            else:
                plt.imshow(image, cmap='gray', vmin=inputvmin, vmax=inputvmax)
            #plt.colorbar()
            plt.gca().invert_yaxis()
    
    plt.subplots_adjust(right=0.8)
    cax = plt.axes([0.85, 0.17, 0.03, 0.7])
    plt.colorbar(cax=cax)
    
#    plt.subplots_adjust(right=0.8)
#    #cbar_ax = plt.add_axes([0.85, 0.15, 0.05, 0.7])
#    plt.colorbar()#cax=cbar_ax)
#    plt.tight_layout()
    
    print 'not_use ',not_use
    #if inputrms == 0 :
    final_image = final_image / ( len(data_tmp.filled())-not_use )    
    return final_image

def PlotImageHist(inputimage,inputN):
    plt.hist(inputimage.flatten(), inputN)
    return
def PlotImage(wcs,inputwcs,inputimage,inputvmin,inputvmax):
    if inputwcs==1:
        plt.subplot(projection=wcs)
    if inputvmax==9999:
        plt.imshow(inputimage, cmap='gray', vmin=inputvmin, vmax=inputimage.max())
    else:
        plt.imshow(inputimage, cmap='gray', vmin=inputvmin, vmax=inputvmax)
    plt.colorbar(label='Jy/beam')
    plt.gca().invert_yaxis()
    return
# =============================================================================
# main code
# =============================================================================

time1 = time.time()



number = 1
catalog = [None]*1
catalog[0] = 'COSMOS2015_Laigle+_v1.1_850wide+850narrow+450narrow+24micron+3GHz_simple.fits' 

###set columns in the main catalog
color1 = "MNUV"
color2 = "MR"
color3 = "MJ"
color = [ color1, color2, color3 ]
photoz = "PHOTOZ"
V_error = "V_MAGERR_APER3"
ip_error = "ip_MAGERR_APER3"
J_error = "J_MAGERR_APER3"
Ks_error = "Ks_MAGERR_APER3"
Ks_mag = "Ks_MAG_APER3"
mass = "MASS_MED"
class_star = "TYPE"
class_SFG = "CLASS"
ra = "ALPHA_J2000"
dec = "DELTA_J2000"


ra = [None]*number
dec = [None]*number
ra[0] = "ALPHA_J2000"
dec[0] = "DELTA_J2000"


###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[0])

######### stack QG without 450 detection

###mask data
mask = []
mask.append((data[0]['450NARROW']==0) & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0) & Mask_myclassQG(0))

###candidate within selected region
#data_tmp = data[0][mask[0]][ MaskRegion('10h00m25.0s', '2d24m22.0s',0.2*u.degree) ]

data_tmp = data[0][mask[0]][ MaskRegion('10h00m30.2s', '02d26m50.0s',7.5*u.arcmin) ]

# stack image #

###read image file
image_data,wcs = ReadImage('04_COSMOS450_850/STUDIES/STUDIES+S2CLS+Casey_450_mf.fits')
#PlotImageHist(image_data[~np.isnan(image_data)],100)
#PlotImage(wcs,1,image_data,-10,10)

###stack
plt.figure(1,figsize=(8,6))
final_image = Stacking(data_tmp,-10,10,0)

plt.figure(2,figsize=(12,5))

plt.subplot(1, 2, 2)
PlotImageHist(final_image,30)
plt.title('450 micron image\n '+r'$\mu$'+str("%.2f" % (np.mean(final_image)))+'Jy/beam, $\sigma$: '\
          +str("%.2f" % (np.std(final_image)))+'Jy/beam')

plt.subplot(1, 2, 1)
PlotImage(wcs,0,final_image,-0.4,0.4)
plt.title('450 micron image')
plt.tight_layout()


'''
# stack rms #

###read image file
image_data,wcs = ReadImage('04_COSMOS450_850/STUDIES/STUDIES+S2CLS+Casey_450_mf_rms.fits')
#PlotImageHist(image_data[~np.isnan(image_data)],100)
#PlotImage(wcs,1,image_data,-10,10)

###stack

plt.figure(4)
final_image = Stacking(0,9999,1)

plt.figure(5)
PlotImageHist(final_image,30)
plt.title('450 micron rms\n mean: '+str("%.2f" % (np.mean(final_image)))+', max: '\
          +str("%.2f" % (np.max(final_image)))+', min: '\
          +str("%.2f" % (np.min(final_image)))+' Jy/beam')

plt.figure(6)
PlotImage(wcs,0,final_image,0,9999)
plt.title('450 micron rms')
'''


######### stack QG with 450 detection

###mask data
mask = []
mask.append((data[0]['450NARROW']!=0) & Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0) & Mask_myclassQG(0))

###candidate within selected region
#data_tmp = data[0][mask[0]][ MaskRegion('10h00m25.0s', '2d24m22.0s',0.2*u.degree) ]

data_tmp = data[0][mask[0]][ MaskRegion('10h00m30.2s', '02d26m50.0s',7.5*u.arcmin) ]


# stack image #

###read image file
image_data,wcs = ReadImage('04_COSMOS450_850/STUDIES/STUDIES+S2CLS+Casey_450_mf.fits')
#PlotImageHist(image_data[~np.isnan(image_data)],100)
#PlotImage(wcs,1,image_data,-10,10)

###stack
plt.figure(3,figsize=(8,6))
final_image = Stacking(data_tmp,-10,10,0)

plt.figure(4,figsize=(12,5))
plt.subplot(1, 2, 2)
PlotImageHist(final_image,30)
plt.title('450 micron image\n '+r'$\mu$'+str("%.2f" % (np.mean(final_image)))+'Jy/beam, $\sigma$: '\
          +str("%.2f" % (np.std(final_image)))+'Jy/beam')

plt.subplot(1, 2, 1)
PlotImage(wcs,0,final_image,-10,10)
plt.title('450 micron image')
plt.tight_layout()


time2 = time.time()
print 'done! time :', time2-time1 , 'sec'