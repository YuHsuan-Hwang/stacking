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
    return (data[index][photoz]>0) & (data[index][photoz]<1.0)

def Mask_class_star(index):
    return (data[index][class_star]==0)| (data[index][class_star]==2)###

def Mask_classQG(index):
    return (data[index][class_SFG]==0)
def Mask_classSFG(index):
    return (data[index][class_SFG]==1)
def Mask_myclassQG(index):
    return (y[index]>3.1) & (y[index]>3.0*x[index]+1.0)
def Mask_myclassSFG(index):
    return (y[index]<=3.1) | (y[index]<=3.0*x[index]+1.0)

def Mask_mass(index,low,high):
    return (data[index][mass]<=high) & (data[index][mass]>low)

def Mask_z(index,low,high):
    return (data[index][photoz]<=high) & (data[index][photoz]>low)

def MeanMass(data):
    return np.mean(10**data[mass])

def MedianRedshift(data):
    return np.ma.median(data[photoz])

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
    if len(data_tmp.filled())>8:
        rnds = random.sample(range(0, len(data_tmp.filled())), 9)
    else:
        rnds = random.sample(range(0, len(data_tmp.filled())), len(data_tmp.filled()))
    print rnds
    cut_num = 0
    for i in range( len(data_tmp.filled()) ):
        ra1 = data_tmp[ra[0]][i]
        dec1 = data_tmp[dec[0]][i]
        x1, y1 = wcs.all_world2pix([[ra1, dec1]], 0)[0]
        imagebox = np.copy(image_data[int(y1)-15:int(y1)+16:1,int(x1)-15:int(x1)+16:1])
        rmsbox = np.copy(rms_data[int(y1)-15:int(y1)+16:1,int(x1)-15:int(x1)+16:1])
        rmssumbox = np.copy(1./rmsbox)
        '''
        for j in range(31):
            for k in range(31):
                if (((imagebox[j][k]/rmsbox[j][k])>=3.5)):
                    imagebox[j][k] = 0
                    rmssumbox[j][k] = 0
                    cut_num +=1
        '''
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
    cbar = plt.colorbar()#label='mJy/beam', fontdict = {'fontsize' : 12})
    cbar.set_label('mJy/beam',size=12)
    plt.gca().invert_yaxis()
    return

def OnePointStacking(data_tmp,wcs,image_data,rms_data):
    middlepix = 0
    middlepix_rmssum = 0
    cut_num = 0
    for i in range( len(data_tmp.filled()) ):
        ra1 = data_tmp[ra[0]][i]
        dec1 = data_tmp[dec[0]][i]
        x1_tmp, y1_tmp = wcs.all_world2pix([[ra1, dec1]], 0)[0]
        x1,y1 = int(x1_tmp),int(y1_tmp)
        
        if (image_data[y1,x1]/rms_data[y1,x1])<SHRESHOLD:
            middlepix = middlepix + image_data[y1,x1]*(1/rms_data[y1,x1])
            middlepix_rmssum = middlepix_rmssum + (1/rms_data[y1,x1])
        else:
            cut_num+=1
        
    middlepix = middlepix/middlepix_rmssum
    return middlepix,cut_num

def PlotRandomPos(inputrandom_ra_list,inputrandom_dec_list):
    plt.scatter(inputrandom_ra_list,inputrandom_dec_list,s=0.1)
    plt.xlabel('ra')
    plt.ylabel('dec')
    plt.axis('scaled')
    return

def task(MASK,TITLE):
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
    
    
    '''
    print 'image stacking'
    
    ###stack
    plt.figure(1,figsize=(8,6))
    final_image = Stacking(data_tmp,wcs,image_data,rms_data,IMAGE_MAX, IMAGE_MIN)
    
    plt.figure(2,figsize=(7,6))
    PlotImage(wcs,0,final_image,STACK_MAX,STACK_MIN)
    plt.title(TITLE+' (N='+str(len(data_tmp.filled()))+')\n'\
              +'mean flux = '+str("%.2f" % (middlepix-MEAN))+'+/-'+str("%.2f" % (SIGMA/np.sqrt(len(data_tmp))))\
              +'mJy/beam, SNR = '+str("%.1f" % ((middlepix-MEAN)/(SIGMA/np.sqrt(len(data_tmp)))))\
              , fontdict = {'fontsize' : 16}, y=1.02)
    plt.xlabel('arcsec', fontdict = {'fontsize' : 12})
    plt.ylabel('arcsec', fontdict = {'fontsize' : 12})
    
    plt.tight_layout()
    
    #fig.savefig('QGstacking.png', format='png', dpi=1200)
    
    print 'Calculation'
    print 'image stacking'
    print 'mean = ',"%.4f" % (np.mean(final_image))
    print 'std = ',"%.4f" % (np.std(final_image))
    print 'SNR : ',"%.4f" % (np.mean(final_image)/np.std(final_image))
    '''
    
    #print 'one point stacking'
    
    middlepix,cut_num = \
    OnePointStacking(data_tmp,wcs,image_data,rms_data)
    
    '''
    print 'middlepix = ',"%.4f" % middlepix
    print 'final result :',"%.4f" % (middlepix-MEAN),'+/-',"%.4f" % (SIGMA/np.sqrt(len(data_tmp)))
    print 'SNR : ',"%.4f" % ((middlepix-MEAN)/SIGMA*np.sqrt(len(data_tmp))),' sigma'
    '''
    mean_flux = middlepix-MEAN
    num = len(data_tmp)-cut_num
    print '{:<18s}{:<9s}{:<10s}{:<6d}{:^8s}{:>6s}'.format(\
           TITLE,("%.1f" % np.log10(MeanMass(data_tmp))),("%.1f" % MedianRedshift(data_tmp)),num,\
           (("%.2f" % mean_flux)+'+/-'+("%.2f" % (SIGMA/np.sqrt(num)))),\
           ("%.1f" % (mean_flux/SIGMA*np.sqrt(num)))  )
    
    return





def RandomStacking(inputpos_num,wcs,image_data,rms_data,inputradius,inputra,inputdec):
    random_ra_list = []
    random_dec_list = []
    test_num = 1000
    pos_num = inputpos_num
    randompix = np.zeros(test_num)
    block_num = 0
    
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
            
            while (1):
                if (image_data[y1,x1]/rms_data[y1,x1])>=SHRESHOLD:

                    
                    block_num = block_num + 1
                    
                    random_ang1 = np.random.uniform(0,2*np.pi,1)
                    random_dis1 = inputradius.to(u.degree).value*np.sqrt(np.random.uniform(0,1,1))
                    random_ra1 = center.ra.value + random_dis*np.cos(random_ang1)
                    random_dec1 = center.dec.value + random_dis*np.sin(random_ang1)
                    
                    x1_tmp, y1_tmp = wcs.all_world2pix([[random_ra1[0], random_dec1[0]]], 0)[0]
                    x1,y1 = int(x1_tmp),int(y1_tmp)

                    
                else:
                    randompix[i] = randompix[i] + image_data[y1,x1]*(1/rms_data[y1,x1])
                    randompix_rmssum = randompix_rmssum + (1/rms_data[y1,x1])
                    break

        randompix[i] = randompix[i]/randompix_rmssum
        
    return randompix, block_num, random_ra_list, random_dec_list

def RandomStackingTask():

    ###read image file
    image_data,wcs = ReadImage(IMAGE,1)
    
    #PlotImageHist(image_data[~np.isnan(image_data)],100)
    #PlotImage(wcs,1,image_data,-10,10)
    
    rms_data = ReadImage(RMS,0)
    
    pos_num_list = [10,100,1000,2000,5000,10000]
    num = len(pos_num_list)
    randompix_mean = np.zeros(num)
    randompix_std = np.zeros(num)
    block_num = np.zeros(num)
    for i in range(num):
        randompix, block_num[i], random_ra_list, random_dec_list = \
        RandomStacking(pos_num_list[i],wcs,image_data,rms_data,RADIUS,CENTER_RA, CENTER_DEC)
        randompix_mean[i] = np.mean(randompix)
        randompix_std[i] = np.std(randompix)#*np.sqrt(pos_num_list[i])
        
        #print 'block_num ',block_num[i]
        #print randompix
        '''
        if i==0:
            plt.figure(1,figsize=(8,6))
            PlotRandomPos(random_ra_list,random_dec_list)
        '''
        #plt.figure(i+1,figsize=(8,6))
        #PlotOnePointStackingHist(0,randompix,30)
        
    
    print randompix_mean 
    plt.figure(1,figsize=(8,6))
    plt.plot(pos_num_list,randompix_mean)
    plt.title('mean')
    plt.show()
    
    print randompix_std
    plt.figure(2,figsize=(8,6))
    plt.plot(pos_num_list,randompix_std)
    plt.plot(pos_num_list,randompix_std[0]*np.sqrt(pos_num_list[0])/np.sqrt(pos_num_list),color="C1")
    #plt.yscale('log')
    plt.title('std')
    plt.show()
    
    sstd = []
    for i in range(num):
        sstd.append(randompix_std[i]*np.sqrt(pos_num_list[i]))
    print sstd
    plt.figure(3,figsize=(8,6))
    plt.plot(pos_num_list,sstd)
    plt.title('sstd')
    plt.show()

def stack_multiband():
    '''
    MASK1 = 
    TITLE = "QGs with 850 $\mu$m detection"
    STACK_MAX,STACK_MIN = -5,5
    
    task(MASK,TITLE)
    MASK = MASK1 & MASK2 & MASK3
    '''
    
    
    MASK1 = (data[0]['450NARROW']==0)
    TITLE = "All QG"
    STACK_MAX,STACK_MIN = -0.1,0.1
    
    MASK = MASK1 & MASK2 & MASK3
    task(MASK,TITLE)
    
    
    
    '''
    MASK1 = (data[0]['450NARROW']==0)&(data[0]['24MICRON']==0)
    TITLE = "QGs without 24 $\mu$m detection"
    STACK_MAX,STACK_MIN = -0.1,0.1
    
    MASK = MASK1 & MASK2 & MASK3
    task(MASK,TITLE)
    '''
    
    
    MASK1 = (data[0]['450NARROW']==0)&(data[0]['24MICRON']==1)
    TITLE = "24-micron ctp."
    STACK_MAX,STACK_MIN = -0.35,0.35
    
    MASK = MASK1 & MASK2 & MASK3
    task(MASK,TITLE)
    
    '''
    MASK1 = (data[0]['450NARROW']==0)&(data[0]['3GHZ']==0)
    TITLE = "QGs without 3 GHz detection"
    STACK_MAX,STACK_MIN = -0.1,0.1
    
    MASK = MASK1 & MASK2 & MASK3
    task(MASK,TITLE)
    '''
    
    
    MASK1 = (data[0]['450NARROW']==0)&(data[0]['3GHZ']==1)
    TITLE = "3-GHz ctp."
    STACK_MAX,STACK_MIN = -0.2,0.2
    
    MASK = MASK1 & MASK2 & MASK3
    task(MASK,TITLE)
    
    
    #MASK1 = (data[0]['450NARROW']==0)&(data[0]['SFG']==1)
    MASK1 = (data[0]['450NARROW']==0)&(data[0]['Radio_excess']==0)
    TITLE = '3-GHz ctp.: SFG'
    STACK_MAX,STACK_MIN = -0.6,0.6
    
    MASK = MASK1 & MASK2 & MASK3
    task(MASK,TITLE)
    
    
    
    #MASK1 = (data[0]['450NARROW']==0)&(data[0]['SFG']==0)
    MASK1 = (data[0]['450NARROW']==0)&(data[0]['Radio_excess']==1)
    TITLE = '3-GHz ctp.: AGN'
    STACK_MAX,STACK_MIN = -0.2,0.2
    
    MASK = MASK1 & MASK2 & MASK3
    task(MASK,TITLE)
    
    
    #MASK1 = (data[0]['450NARROW']==0)&(data[0]['24MICRON']==0)&(data[0]['SFG']!=1)
    MASK1 = (data[0]['450NARROW']==0)&(data[0]['24MICRON']==0)&(data[0]['Radio_excess']!=0)
    TITLE = "'QG'"
    STACK_MAX,STACK_MIN = -0.1,0.1
    
    MASK = MASK1 & MASK2 & MASK3
    task(MASK,TITLE)
    
    
    #MASK1 = (data[0]['450NARROW']==0)&( (data[0]['24MICRON']==1)|(data[0]['SFG']==1) )
    MASK1 = (data[0]['450NARROW']==0)&( (data[0]['24MICRON']==1)|(data[0]['Radio_excess']==0) )
    TITLE = "IR-bright QG'"
    STACK_MAX,STACK_MIN = -0.1,0.1
    
    MASK = MASK1 & MASK2 & MASK3
    task(MASK,TITLE)
    
    
    

    
    return

def stack_mass(mass_cut):
    
    MASK1 = (data[0]['450NARROW']==0)&(data[0]['24MICRON']==0)&(data[0]['Radio_excess']!=0)
    
    #QGs without "dusty SFGs" small mass
    TITLE = 'logM* < '+str(mass_cut)
    MASK_mass = Mask_mass(0,0,mass_cut)
    STACK_MAX,STACK_MIN = -0.1,0.1
    MASK = MASK_mass & MASK1 & MASK2 & MASK3
    task(MASK,TITLE)
    
    #QGs without "dusty SFGs" large mass
    TITLE = 'logM* > '+str(mass_cut)
    MASK_mass = Mask_mass(0,mass_cut,100)
    STACK_MAX,STACK_MIN = -0.1,0.1
    MASK = MASK_mass & MASK1 & MASK2 & MASK3
    task(MASK,TITLE)
    return

def stack_z(z_cut):
    
    MASK1 = (data[0]['450NARROW']==0)&(data[0]['24MICRON']==0)&(data[0]['Radio_excess']!=0)
    
    #QGs without "dusty SFGs" small mass
    TITLE = 'z < '+str(z_cut)
    MASK_z = Mask_z(0,0,z_cut)
    STACK_MAX,STACK_MIN = -0.1,0.1
    MASK = MASK_z & MASK1 & MASK2 & MASK3
    task(MASK,TITLE)
    
    #QGs without "dusty SFGs" large mass
    TITLE = 'z > '+str(z_cut)
    MASK_z = Mask_z(0,z_cut,100)
    STACK_MAX,STACK_MIN = -0.1,0.1
    MASK = MASK_z & MASK1 & MASK2 & MASK3
    task(MASK,TITLE)
    return

def stackPrintBin():
    
    stack_nondustySFG()
    print
    
    #mass
    mass_cut = 10
    stack_mass(mass_cut)
    print
    mass_cut = 10.5
    stack_mass(mass_cut)
    print
    mass_cut = 11
    stack_mass(mass_cut)
    print
    mass_cut = 11.5
    stack_mass(mass_cut)
    '''
    '''
    #redshift
    z_cut = 0.5
    stack_z(z_cut)
    print
    z_cut = 0.75
    stack_z(z_cut)
    print
    z_cut = 1
    stack_z(z_cut)
    print
    z_cut = 1.5
    stack_z(z_cut)
    print
    z_cut = 2
    stack_z(z_cut)

# =============================================================================
# main code
# =============================================================================

time1 = time.time()

number = 1
catalog = [None]*1
catalog[0] = "COSMOS2015_merged.fits"

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

CENTER_RA,CENTER_DEC = '10h00m25.0s', '2d24m22.0s'
#'10h00m30.2s', '02d26m50.0s' #'10h00m25.0s', '2d24m22.0s'
RADIUS = 0.2*u.degree #7.5*u.arcmin #0.2*u.degree

path = "/Users/yuhsuan/Documents/research/05WH/data/COSMOS/"

###read catalog
data = [None]*number
x = [None]*number
y = [None]*number

for i in range(number):
    ReadCatalog(i,catalog[0])

###set shreshold
SHRESHOLD = 3.0

######### random stacking
'''
print 'random stacking'

MASK2 = Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0)
MASK3 = Mask_myclassQG(0)
MASK = MASK2 & MASK3

IMAGE = path+'04_COSMOS450_850/STUDIES/STUDIES+S2CLS+Casey_450_mf.fits'
RMS = path+'04_COSMOS450_850/STUDIES/STUDIES+S2CLS+Casey_450_mf_rms.fits'

RandomStackingTask()
'''

######### stack QG


if (SHRESHOLD==10000.0):
    MEAN = 0.000  #0.0
    SIGMA = 2.414  #2.4
elif (SHRESHOLD==3.5):
    MEAN = -0.126
    SIGMA = 2.073
elif (SHRESHOLD==3.0):
    MEAN = -0.153
    SIGMA = 2.005
elif (SHRESHOLD==2.0):
    MEAN = -0.290
    SIGMA = 1.931
if (SHRESHOLD==10000.0):
    print "shreshold = -    mean =",MEAN,"   sigma =",SIGMA,"/sqrt(N)"
else:
    print "shreshold =",SHRESHOLD,"    mean =",MEAN,"   sigma =",SIGMA,"/sqrt(N)"


#print "450 micron stacking"

print "-"*85
#print '{:<30s}{:<20s}{:^12s}{:>10s}'.format("Groups","Number","Mean flux","SNR" )
print '{:<16s}{:<10s}{:<10s}{:<8s}{:^10s}{:>6s}'.format("Groups","Mean_M*","Median_z","Number","Mean flux","SNR" )

print "-"*85

MASK2 = Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0)
MASK3 = Mask_myclassQG(0)

IMAGE_MAX, IMAGE_MIN = -10,10

IMAGE = path+'04_COSMOS450_850/STUDIES/STUDIES+S2CLS+Casey_450_mf.fits'
RMS = path+'04_COSMOS450_850/STUDIES/STUDIES+S2CLS+Casey_450_mf_rms.fits'

stack_multiband()

#stackPrintBin()


print "-"*85
print



time2 = time.time()
print 'done! time :', time2-time1 , 'sec'
