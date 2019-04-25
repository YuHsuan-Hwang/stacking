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

def MaskRegion(mask,inputra,inputdec,inputradius):
    c0 = SkyCoord(ra=data[0][mask[0]][ra[0]], dec=data[0][mask[0]][dec[0]])
    center = SkyCoord(inputra, inputdec, frame='fk5')
    radius = inputradius
    sep = center.separation(c0)
    return sep<=radius

def ReadImage(filename,wcs_flag):
    hdu_list = fits.open(filename)
    #hdu_list.info()
    data = hdu_list[0].data[0]
    if wcs_flag==1:
        wcs = WCS(hdu_list[0].header).celestial
        return data,wcs
    else:
        return data

def Stacking(data_tmp,wcs,image_data,rms_data,inputvmin,inputvmax):
    
    final_image = np.zeros((31,31))
    rms_sum = np.zeros((31,31))
    not_use = 0
    print 'len(data_tmp.filled()) ',len(data_tmp.filled())
   
    picnum = 1
    
    rnds = random.sample(range(0, len(data_tmp.filled())), 9) 
    print rnds
    cut_num = 0
    for i in range( len(data_tmp.filled()) ):
        ra1 = data_tmp[ra[0]][i]
        dec1 = data_tmp[dec[0]][i]
        x1, y1 = wcs.all_world2pix([[ra1, dec1]], 0)[0]
        imagebox = np.copy(image_data[int(y1)-15:int(y1)+16:1,int(x1)-15:int(x1)+16:1])
        rmsbox = np.copy(rms_data[int(y1)-15:int(y1)+16:1,int(x1)-15:int(x1)+16:1])
        rmssumbox = np.copy(1./rmsbox)

        for j in range(31):
            for k in range(31):
                if (((imagebox[j][k]/rmsbox[j][k])>=3.5)):
                #& ( (j<11)|(j>18)|(k<11)|(k>18) )
                    #& ( (j<13)|(j>16)|(k<13)|(k>16) )
                    #|(imagebox[j][k]/rmsbox[j][k]<=-3.5)
                    imagebox[j][k] = 0
                    rmssumbox[j][k] = 0
                    cut_num +=1

        
        image = imagebox*rmssumbox
        rms_sum += rmssumbox
        
        if (image.shape==(31,31)) & (np.isnan(image).any()==False):
            final_image += image
        else:
            not_use += 1
            
        if i in rnds:
            plt.subplot(3, 3, picnum)
            picnum = picnum+1
            if inputvmax==9999:
                plt.imshow(imagebox, cmap='gray', vmin=inputvmin, vmax=image.max())
            else:
                plt.imshow(imagebox, cmap='gray', vmin=inputvmin, vmax=inputvmax)
            #plt.colorbar()
            plt.gca().invert_yaxis()
            plt.title(i)
    plt.tight_layout()
    plt.subplots_adjust(right=0.8)
    cax = plt.axes([0.85, 0.17, 0.03, 0.7])
    plt.colorbar(cax=cax)
    
    print 'not_use ',not_use
    print 'cut_num ',cut_num
    #final_image = final_image / ( len(data_tmp.filled())-not_use )
    final_image = final_image / rms_sum
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

def OnePointStacking(data_tmp,wcs,image_data,rms_data,inputradius,inputra,inputdec):
    
    middlepix = 0
    middlepix_rmssum = 0
    for i in range( len(data_tmp.filled()) ):
        ra1 = data_tmp[ra[0]][i]
        dec1 = data_tmp[dec[0]][i]
        x1_tmp, y1_tmp = wcs.all_world2pix([[ra1, dec1]], 0)[0]
        x1,y1 = int(x1_tmp),int(y1_tmp)
        if (image_data[y1,x1]/rms_data[y1,x1])<3.5:
            middlepix = middlepix + image_data[y1,x1]*(1/rms_data[y1,x1])
            middlepix_rmssum = middlepix_rmssum + (1/rms_data[y1,x1])
    middlepix = middlepix/middlepix_rmssum
    
    random_ra_list = []
    random_dec_list = []
    test_num = 1000
    pos_num = len(data_tmp.filled())
    randompix = np.zeros(test_num)
    for i in range(test_num):
        random_ang = np.random.uniform(0,2*np.pi,pos_num)
        random_dis = inputradius.to(u.degree).value*np.sqrt(np.random.uniform(0,1,pos_num))
        center = SkyCoord(inputra, inputdec, frame='fk5')
        random_ra = center.ra.value + random_dis*np.cos(random_ang)
        random_dec = center.dec.value + random_dis*np.sin(random_ang)
        if i == 0:
            random_ra_list.append(random_ra)
            random_dec_list.append(random_dec)
        randompix_rmssum = 0
        for j in range(pos_num):
            ra1 = random_ra[j]
            dec1 = random_dec[j]
            x1_tmp, y1_tmp = wcs.all_world2pix([[ra1, dec1]], 0)[0]
            x1,y1 = int(x1_tmp),int(y1_tmp)
            if (image_data[y1,x1]/rms_data[y1,x1])<3.5:
                randompix[i] = randompix[i] + image_data[y1,x1]*(1/rms_data[y1,x1])
                randompix_rmssum = randompix_rmssum + (1/rms_data[y1,x1])
        randompix[i] = randompix[i]/randompix_rmssum
        
    return middlepix, randompix, random_ra_list, random_dec_list

def PlotOnePointStackingHist(inputmiddlepix,inputrandompix,inputN):
    plt.hist(inputrandompix, inputN)
    plt.title(IMAGE_NAME+': '+r'$\mu$:'+str("%.4f" % (inputmiddlepix))+'\n'\
              +r'$\mu$:'+str("%.4f" % (np.mean(inputrandompix)))+'mJy/beam, $\sigma$: '\
              +str("%.4f" % (np.std(inputrandompix)))+'mJy/beam')

def PlotRandomPos(inputrandom_ra_list,inputrandom_dec_list):
    plt.scatter(inputrandom_ra_list,inputrandom_dec_list,s=0.1)
    plt.xlabel('ra')
    plt.ylabel('dec')
    plt.axis('scaled')
    return

def task():
    ###mask data
    mask = []
    mask.append(MASK)
    
    ###candidate within selected region
    data_tmp = data[0][mask[0]][ MaskRegion(mask,CENTER_RA, CENTER_DEC, RADIUS) ]
    
    # stack image #
    
    ###read image file
    image_data,wcs = ReadImage(IMAGE,1)
    #PlotImageHist(image_data[~np.isnan(image_data)],100)
    #PlotImage(wcs,1,image_data,-10,10)
    rms_data = ReadImage(RMS,0)
    
    ###stack
    plt.figure(1,figsize=(8,6))
    final_image = Stacking(data_tmp,wcs,image_data,rms_data,IMAGE_MAX, IMAGE_MIN)
    
    plt.figure(2,figsize=(12,5))
    
    plt.subplot(1, 2, 2)
    PlotImageHist(final_image,30)
    plt.title(IMAGE_NAME+'\n '+r'$\mu$:'+str("%.4f" % (np.mean(final_image)))+'mJy/beam, $\sigma$: '\
              +str("%.4f" % (np.std(final_image)))+'mJy/beam')
    
    plt.subplot(1, 2, 1)
    PlotImage(wcs,0,final_image,STACK_MAX,STACK_MIN)
    plt.title(IMAGE_NAME+'\n '+str(len(data_tmp.filled()))+' galaxies')
    plt.tight_layout()
    
    middlepix, randompix,random_ra_list,random_dec_list = \
    OnePointStacking(data_tmp,wcs,image_data,rms_data,RADIUS,CENTER_RA, CENTER_DEC)
    
    plt.figure(3,figsize=(8,6))
    PlotOnePointStackingHist(middlepix,randompix,30)
    
    plt.figure(4,figsize=(8,6))
    PlotRandomPos(random_ra_list,random_dec_list)
    
    print 'Calculation'
    print 'middlepix = ',"%.4f" % middlepix
    print 'randompix_mean = ',"%.4f" % np.mean(randompix)
    print 'randompix_std = ',"%.4f" % np.std(randompix)
    print 'final result :',"%.4f" % (middlepix-np.mean(randompix)),'+/-'\
    ,"%.4f" % np.std(randompix)
    print 'SNR : ',"%.4f" % ((middlepix-np.mean(randompix))/np.std(randompix)),' sigma'
    
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

print 'stack QG without 450 detection'
MASK1 = (data[0]['450NARROW']==0) & (data[0]['24MICRON']==1) & (data[0]['3GHZ']==0)
# (data[0]['24MICRON']==1)#& (data[0]['3GHZ']==1)
MASK2 = Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0)
MASK3 = Mask_myclassQG(0)
MASK = MASK1 & MASK2 & MASK3

CENTER_RA,CENTER_DEC = '10h00m25.0s', '2d24m22.0s'
#'10h00m30.2s', '02d26m50.0s' #'10h00m25.0s', '2d24m22.0s'
RADIUS = 0.2*u.degree #7.5*u.arcmin #0.2*u.degree

IMAGE_MAX, IMAGE_MIN = -10,10
STACK_MAX,STACK_MIN = -1,1#-0.4,0.4#


IMAGE = '04_COSMOS450_850/STUDIES/STUDIES+S2CLS+Casey_450_mf.fits'
IMAGE_NAME = '450 micron image'

RMS = '04_COSMOS450_850/STUDIES/STUDIES+S2CLS+Casey_450_mf_rms.fits'

task()



######### stack QG with 450 detection
'''
print 'stack QG with 450 detection'
MASK1 = data[0]['450NARROW']!=0
MASK2 = Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0)
MASK3 = Mask_myclassQG(0)
MASK = MASK1 & MASK2 & MASK3

CENTER_RA,CENTER_DEC = '10h00m25.0s', '2d24m22.0s'
#'10h00m30.2s', '02d26m50.0s' #'10h00m25.0s', '2d24m22.0s'
RADIUS = 0.2*u.degree #7.5*u.arcmin #0.2*u.degree

IMAGE_MAX, IMAGE_MIN = -10,10
STACK_MAX,STACK_MIN = -10,10


IMAGE = '04_COSMOS450_850/STUDIES/STUDIES+S2CLS+Casey_450_mf.fits'
IMAGE_NAME = '450 micron image'

RMS = '04_COSMOS450_850/STUDIES/STUDIES+S2CLS+Casey_450_mf_rms.fits'

task()
'''

time2 = time.time()
print 'done! time :', time2-time1 , 'sec'