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
from astropy.cosmology import FlatLambdaCDM
# =============================================================================
# functions
# =============================================================================

def ReadCatalog( index, catalog ):
    
    data[index] = Table.read( catalog, hdu=1 )
    x.append( data[index][color2]-data[index][color3] )
    y.append( data[index][color1]-data[index][color2] )
    
    return


def Mask_M(index):
    
    down, up    = -30, 0
    mask_color1 = ( data[index][color1]>down ) & ( data[index][color1]<up )
    mask_color2 = ( data[index][color2]>down ) & ( data[index][color2]<up )
    mask_color3 = ( data[index][color3]>down ) & ( data[index][color3]<up )
   
    return mask_color1 & mask_color2 & mask_color3

def Mask_error( mode, up, index ):
    
    if mode==1:
        mask_V_errorV = ( data[index][V_error] <up ) & ( data[index][V_error] >0 )
        mask_ip_error = ( data[index][ip_error]<up ) & ( data[index][ip_error]>0 ) #Suprime-Cam:i+band
        mask_J_error  = ( data[index][J_error] <up ) & ( data[index][J_error] >0 )
        mask_Ks_error = ( data[index][Ks_error]<up ) & ( data[index][Ks_error]>0 )
        return mask_V_errorV | mask_ip_error | mask_J_error | mask_Ks_error
    
    if mode==2:
        return (data[index][Ks_error]<up) & (data[index][Ks_error]>0)
    
    if mode==3:
        return (data[index][Ks_mag]<24)
    

def Mask_photoz( index ):
    return (data[index][photoz]>0)

def Mask_class_star( index ):
    return ( data[index][class_star]==0 ) | ( data[index][class_star]==2 )

def MaskAll( index ):
    return Mask_M( index ) & Mask_error( 1, 0.1, index ) & Mask_photoz( index ) & Mask_class_star( index )
def MaskAll_oneband( index ):
    return Mask_M( index ) & Mask_error( 2, 0.2, index ) & Mask_photoz( index ) & Mask_class_star( index )



def Mask_classQG   ( index ): return ( data[index][class_SFG]==0 )
def Mask_classSFG  ( index ): return ( data[index][class_SFG]==1 )
def Mask_myclassQG ( index ): return ( y[index]> 3.1 ) & ( y[index]> 3.0*x[index]+1.0 )
def Mask_myclassSFG( index ): return ( y[index]<=3.1 ) | ( y[index]<=3.0*x[index]+1.0 )

def MaskAllQG  ( index ): return MaskAll( index ) & Mask_myclassQG( index )
def MaskAllQG_oneband  ( index ): return MaskAll_oneband( index ) & Mask_myclassQG( index )

def Mask_nonSMG_850( index ): return ( data[index]['850WIDE'  ] ==0 ) | ( data[index]['850WIDE'  ] ==7 )
def Mask_nonSMG_450( index ): return ( data[index]['450NARROW'] ==0 ) | ( data[index]['450NARROW'] ==7 )
def Mask_nonSMG_850_oneband( index ): return ( data[index]['850WIDE_ONEBAND'  ] ==0 ) | ( data[index]['850WIDE_ONEBAND'  ] ==7 )
def Mask_nonSMG_450_oneband( index ): return ( data[index]['450NARROW_ONEBAND'] ==0 ) | ( data[index]['450NARROW_ONEBAND'] ==7 )


def Mask_nonIRB( index ): return ( data[index]['24MICRON']==0 ) & ( data[index]['Radio_excess']!=0 ) #radio excess = -99 or 1
def Mask_IRB   ( index ): return ( data[index]['24MICRON']==1 ) | ( data[index]['Radio_excess']==0 ) 

def Mask_mass( index, low, high ): return ( data[index][mass  ]<=high ) & ( data[index][mass  ]>low )
def Mask_z   ( index, low, high ): return ( data[index][photoz]<=high ) & ( data[index][photoz]>low )

def MaskRegion( mask,inputra,inputdec,inputradius ): #mask is not list
    
    c0     = SkyCoord( ra=data[0][mask][ra[0]], dec=data[0][mask][dec[0]] )
    center = SkyCoord( inputra, inputdec, frame='fk5' )
    radius = inputradius
    sep    = center.separation( c0 )
    
    return sep<=radius

def MaskRegion_nomask( inputra,inputdec,inputradius ): 
    
    c0     = SkyCoord( ra=data[0][ra[0]], dec=data[0][dec[0]] )
    center = SkyCoord( inputra, inputdec, frame='fk5' )
    radius = inputradius
    sep    = center.separation( c0 )
    
    return sep<=radius

def MeanMass      ( data ):  return np.mean     ( 10**data[mass  ] )
def MedianRedshift( data ):  return np.ma.median(     data[photoz] )

def ReadImage( filename, wcs_flag ):
    
    hdu_list = fits.open(filename)
    #hdu_list.info()
    data = hdu_list[0].data[0]
    if ( wcs_flag==1 ):
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



def OnePointStacking( data_tmp, wcs, image_data, rms_data ):
    
    middlepix        = 0
    middlepix_rmssum = 0
    cut_num          = 0
    
    for i in range( len(data_tmp.filled()) ):
        
        ra1            = data_tmp[ra[0] ][i]
        dec1           = data_tmp[dec[0]][i]
        x1_tmp, y1_tmp = wcs.all_world2pix( [[ra1, dec1]], 0 )[0]
        x1    , y1     = int( x1_tmp ), int( y1_tmp )
        
        if ( (image_data[y1,x1]/rms_data[y1,x1])<SHRESHOLD ):
            middlepix = middlepix + image_data[y1,x1]*(1/rms_data[y1,x1])
            middlepix_rmssum = middlepix_rmssum + (1/rms_data[y1,x1])
            
        else: cut_num+=1
        
    middlepix = middlepix/middlepix_rmssum
    
    return middlepix, cut_num



def PlotRandomPos(inputrandom_ra_list,inputrandom_dec_list):
    plt.scatter(inputrandom_ra_list,inputrandom_dec_list,s=0.1)
    plt.xlabel('ra')
    plt.ylabel('dec')
    plt.axis('scaled')
    return

def task(TITLE,MASK):
    ###mask data
    mask = []
    mask.append(MASK)
    
    ###candidate within selected region
    data_tmp = data[0][mask[0]]
    
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
    
    print 'mean = ',"%.4f" % (np.mean(final_image))
    print 'std = ',"%.4f" % (np.std(final_image))
    print 'SNR : ',"%.4f" % (np.mean(final_image)/np.std(final_image))
    '''
    
    #print 'one point stacking'
    
    middlepix = \
    OnePointStacking(data_tmp,wcs,image_data,rms_data)
    
    print 'middlepix = ',"%.4f" % middlepix
    print 'final result :',"%.4f" % (middlepix-MEAN),'+/-',"%.4f" % (SIGMA/np.sqrt(len(data_tmp)))
    print 'SNR : ',"%.4f" % ((middlepix-MEAN)/SIGMA*np.sqrt(len(data_tmp))),' sigma'
    
    return

def taskprint(TITLE,MASK):
    
    ###candidate within selected region
    data_tmp = data[0][MASK]  
    
    ###stack image  
    middlepix, cut_num = OnePointStacking( data_tmp, wcs, image_data, rms_data )
    
    mean_flux = middlepix - MEAN
    num       = len( data_tmp ) - cut_num
    median_z  = MedianRedshift( data_tmp )
    mean_mass = MeanMass      ( data_tmp )
    
    '''
    mean_mass = MeanMass(data_tmp)
    global LIR_iter
    L_IR_fl = LIR_WH[LIR_iter]
    if (LIR_limitmask_WH[LIR_iter]==1):
        L_IR = "<"+("%.1E" %  L_IR_fl)
    else:
        L_IR = ("%.1E" % L_IR_fl)
    
    LIR_iter+=1
    
    SFR = L_IR_fl*1.7e-10
    sSFR = np.log10(SFR/mean_mass)  
    '''
    
    if (mean_flux<=0.0):
       #### mean_mass = 0.0 
        L_IR = 0.0
        SFR = 0.0
        sSFR = 0.0
    else:
        L_IR = mean_flux*1.9e12
        SFR  = L_IR*1.7e-10
        sSFR = np.log10( SFR/mean_mass )
    
    
    print '{:<18s}{:<9s}{:<10s}{:<6d}{:^8s}{:>6s}{:>12s}{:>6s}{:>6s}'.format(
           TITLE,
           ( "%.1f" % np.log10(MeanMass(data_tmp))   ),
           ( "%.1f" % median_z                       ),
           num,
           ( ("%.2f" % mean_flux) + '+/-' + ("%.2f" % (SIGMA/np.sqrt(num))) ),
           ( "%.1f" % (mean_flux/SIGMA*np.sqrt(num)) ),
           ( "%.1e" % L_IR                           ),
           ( "%3d"  % SFR                            ),
           ( "%.1f" % sSFR                           )  )          

    
    return

def taskprint_latex(TITLE,MASK):
    
    ###candidate within selected region
    data_tmp = data[0][MASK]  
    
    ###stack image  
    middlepix, cut_num = OnePointStacking( data_tmp, wcs, image_data, rms_data )
    
    mean_flux = middlepix - MEAN
    num       = len( data_tmp ) - cut_num
    median_z  = MedianRedshift( data_tmp )
    mean_mass = MeanMass      ( data_tmp )
    
    global LIR_iter
        
    L_IR = LIR_WH[LIR_iter]
    
    LIR_iter+=1
    
    SFR = L_IR*1.7e-10
    
    if (mean_flux<=0.0):
       #### mean_mass = 0.0 
        L_IR = 0.0
        SFR = 0.0
    
    print TITLE \
          +" & "+( "%.1f" % np.log10( mean_mass )   ) \
          +" & "+( "%.1f" % median_z   ) \
          +" & "+str( num ) \
          +" & "+("%.2f" % mean_flux) + "$\pm$" + ("%.2f" % (SIGMA/np.sqrt(num)))  \
          +" & "+( "%.1f" % (mean_flux/SIGMA*np.sqrt(num)) )\
          +" & "+("%.1E" %  L_IR) \
          +" & "+( "%.1f"  % SFR  )+"\\\\"         

    #All QG                 & 10.7 & 0.9 & 17535 & 0.06$\pm$0.01 & 6.2 & 1.0$\times10^{10}$ & 1.8\\
    return
    

def taskSSFR(MASK):
    
    ###candidate within selected region
    data_tmp = data[0][MASK]
    
    ###stack image
    middlepix = OnePointStacking( data_tmp, wcs, image_data, rms_data )
    
    mean_flux = middlepix - MEAN
    
    if (mean_flux<=0.0):
        mean_mass = 0.0
        SFR = 0.0
    else:
        mean_mass = np.log10( MeanMass(data_tmp) )
        SFR = np.log10( (mean_flux*1.9*1.7*100) )

    return mean_mass, SFR

def taskSSFR_WH(MASK):
    
    data_tmp = data[0][MASK]
    
    global LIR_iter
    
    mean_mass = np.log10( MeanMass(data_tmp) )
    SFR = np.log10( LIR_WH[LIR_iter]*1.7e-10 )
    SFRerr_up = np.log10( (LIR_WH[LIR_iter]+LIR_errlen_WH[LIR_iter])*1.7e-10 )\
                - np.log10( LIR_WH[LIR_iter]*1.7e-10 )
    if (LIR_WH[LIR_iter]-LIR_errlen_WH[LIR_iter])<=0:
        SFRerr_low = 100
    else:
        SFRerr_low = np.log10( LIR_WH[LIR_iter]*1.7e-10 )\
                     - np.log10( (LIR_WH[LIR_iter]-LIR_errlen_WH[LIR_iter])*1.7e-10 )
    print LIR_iter, LIR_WH[LIR_iter], SFR
    LIR_iter+=1
    
    return mean_mass, SFR, SFRerr_low, SFRerr_up

def RandomStacking(inputpos_num,wcs,image_data,rms_data,ra_min,ra_max,dec_min,dec_max):
    random_ra_list = []
    random_dec_list = []
    test_num = 1000
    pos_num = inputpos_num
    randompix = np.zeros(test_num)
    block_num = 0
    
    for i in range(test_num):
        
        random1 = np.random.uniform(0,1,pos_num)
        random2 = np.random.uniform(0,1,pos_num)
        random_ra = random1*(ra_max-ra_min) + ra_min
        random_dec = random2*(dec_max-dec_min) + dec_min
        
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
                if ( (image_data[y1,x1]/rms_data[y1,x1])>=SHRESHOLD )|(np.isnan(image_data[y1,x1])):

                    
                    block_num = block_num + 1
                    
                    onerandom1 = np.random.uniform(0,1,1)
                    onerandom2 = np.random.uniform(0,1,1)
                    onerandom_ra = onerandom1*(ra_max-ra_min) + ra_min
                    onerandom_dec = onerandom2*(dec_max-dec_min) + dec_min
                    
                    x1_tmp, y1_tmp = wcs.all_world2pix([[onerandom_ra[0], onerandom_dec[0]]], 0)[0]
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
    
    mask = []
    mask.append(MASK)
    ra_min = data[0][mask[0]][ra[0]].min()
    ra_max = data[0][mask[0]][ra[0]].max()
    dec_min = data[0][mask[0]][dec[0]].min()
    dec_max = data[0][mask[0]][dec[0]].max()
    #print ra_min, ra_max, dec_min, dec_max

    pos_num_list = [10,100,1000,2000,5000,10000]
    num = len(pos_num_list)
    randompix_mean = np.zeros(num)
    randompix_std = np.zeros(num)
    block_num = np.zeros(num)
    for i in range(num):
        randompix, block_num[i], random_ra_list, random_dec_list = \
        RandomStacking(pos_num_list[i],wcs,image_data,rms_data,ra_min,ra_max,dec_min,dec_max)
        randompix_mean[i] = np.mean(randompix)
        randompix_std[i] = np.std(randompix)#*np.sqrt(pos_num_list[i])
        
        #print 'block_num ',block_num[i]
        #print randompix
        '''
        if i==1:
            plt.figure(1,figsize=(8,6))
            PlotRandomPos(random_ra_list,random_dec_list)
        '''
        '''
        plt.figure(i+1,figsize=(8,6))
        plt.hist(randompix,30)
        plt.show
        '''
    print randompix_mean 
    plt.figure(num+1,figsize=(8,6))
    plt.plot(pos_num_list,randompix_mean)
    plt.title('mean')
    plt.show()
    
    print randompix_std
    plt.figure(num+2,figsize=(8,6))
    plt.plot(pos_num_list,randompix_std)
    plt.plot(pos_num_list,randompix_std[0]*np.sqrt(pos_num_list[0])/np.sqrt(pos_num_list),color="C1")
    #plt.yscale('log')
    plt.title('std')
    plt.show()
    
    sstd = []
    for i in range(num):
        sstd.append(randompix_std[i]*np.sqrt(pos_num_list[i]))
    print sstd
    plt.figure(num+3,figsize=(8,6))
    plt.plot(pos_num_list,sstd)
    plt.title('sstd')
    plt.show()



def StackBin( mask ):
    
    #z_list = [ "z <= 0.5", "0.5 < z <= 1.0", "1.0 < z <= 1.5", "1.5 < z <= 2.0", "z > 2.0" ]
    #z_mask_list = [ Mask_z(0,0,0.5), Mask_z(0,0.5,1.0), Mask_z(0,1.0,1.5), 
    #                Mask_z(0,1.5,2.0), Mask_z(0,2.0,100.0) ]
    z_list = [ "z <= 0.5", "0.5 < z <= 1.0", "1.0 < z <= 1.5", "1.5 < z <= 2.0", "2.0 < z <= 2.5", "z > 2.5" ]
    z_mask_list = [ Mask_z(0,0,0.5), Mask_z(0,0.5,1.0), Mask_z(0,1.0,1.5), 
                    Mask_z(0,1.5,2.0), Mask_z(0,2.0,2.5), Mask_z(0,2.5,100.0) ]
    for i in range( len(z_list) ):
        print
        print z_list[i]
        taskprint_latex( 'logM* <= 10.5', mask & z_mask_list[i] & Mask_mass(0,0   ,10.5) )
        taskprint_latex( 'logM* > 10.5',  mask & z_mask_list[i] & Mask_mass(0,10.5,100 ) )
    
    return

def PrintStacking():
    
    if (SHRESHOLD==10000.0):
        print "shreshold = -    mean =",MEAN,"   sigma =",SIGMA,"/sqrt(N)"
    else:
        print "shreshold =",SHRESHOLD,"    mean =",MEAN,"   sigma =",SIGMA,"/sqrt(N)"
    
    dashnum = 85
    
    print "-"*dashnum
    print '{:<16s}{:<10s}{:<10s}{:<8s}{:^10s}{:>6s}{:>10s}{:>8s}{:>6s}'.format(
    "Groups","Mean_M*","Median_z","Number","Mean_flux","SNR","L_IR","SFR","sSFR")
    
    print "-"*dashnum
    
    
    #QGs without 850um detection    #STACK_MAX,STACK_MIN = -0.1,0.1
    taskprint( "All QG",          MaskAllQG(0) & Mask_nonSMG_850(0)                                )
    
    #QGs with 24um detection        #STACK_MAX,STACK_MIN = -0.35,0.35
    taskprint( "24-micron ctp.",  MaskAllQG(0) & Mask_nonSMG_850(0) & (data[0]['24MICRON']    ==1) )
    
    #QGs with 3GHz detection        #STACK_MAX,STACK_MIN = -0.2,0.2
    taskprint( "3-GHz ctp.",      MaskAllQG(0) & Mask_nonSMG_850(0) & (data[0]['3GHZ']        ==1) )
    
    #QGs with label "SFG"           #STACK_MAX,STACK_MIN = -0.6,0.6
    taskprint( "3-GHz ctp.: SFG", MaskAllQG(0) & Mask_nonSMG_850(0) & (data[0]['Radio_excess']==0) )
    
    #QGs with label "AGN"           #STACK_MAX,STACK_MIN = -0.2,0.2
    taskprint( "3-GHz ctp.: AGN", MaskAllQG(0) & Mask_nonSMG_850(0) & (data[0]['Radio_excess']==1) )  
    
    
    print "-"*dashnum
    
    maskQG = MaskAllQG(0) & Mask_nonSMG_850(0) & Mask_nonIRB(0)
    taskprint( "IR-faint QG", maskQG )
    StackBin( maskQG )  

    
    print "-"*dashnum
    
    maskIRBQG = MaskAllQG(0) & Mask_nonSMG_850(0) & Mask_IRB(0)
    taskprint( "IR-bright QG", maskIRBQG )
    StackBin( maskIRBQG )
    
    
    print "-"*dashnum
    print
    
    return

def PrintStackingOneBand():
    
    if (SHRESHOLD==10000.0):
        print "shreshold = -    mean =",MEAN,"   sigma =",SIGMA,"/sqrt(N)"
    else:
        print "shreshold =",SHRESHOLD,"    mean =",MEAN,"   sigma =",SIGMA,"/sqrt(N)"
    
    dashnum = 85
    
    print "-"*dashnum
    print '{:<16s}{:<10s}{:<10s}{:<8s}{:^10s}{:>6s}{:>10s}{:>8s}{:>6s}'.format(
    "Groups","Mean_M*","Median_z","Number","Mean_flux","SNR","L_IR","SFR","sSFR")
    
    print "-"*dashnum
    
    
    #QGs without 850um detection    #STACK_MAX,STACK_MIN = -0.1,0.1
    taskprint_latex( "All QG",          MaskAllQG_oneband(0) & Mask_nonSMG_850_oneband(0)                                )
    
    #QGs with 24um detection        #STACK_MAX,STACK_MIN = -0.35,0.35
    taskprint_latex( "24-micron ctp.",  MaskAllQG_oneband(0) & Mask_nonSMG_850_oneband(0) & (data[0]['24MICRON']    ==1) )
    
    #QGs with 3GHz detection        #STACK_MAX,STACK_MIN = -0.2,0.2
    taskprint_latex( "3-GHz ctp.",      MaskAllQG_oneband(0) & Mask_nonSMG_850_oneband(0) & (data[0]['3GHZ']        ==1) )
    
    #QGs with label "SFG"           #STACK_MAX,STACK_MIN = -0.6,0.6
    taskprint_latex( "3-GHz ctp.: SFG", MaskAllQG_oneband(0) & Mask_nonSMG_850_oneband(0) & (data[0]['Radio_excess']==0) )
    
    #QGs with label "AGN"           #STACK_MAX,STACK_MIN = -0.2,0.2
    taskprint_latex( "3-GHz ctp.: AGN", MaskAllQG_oneband(0) & Mask_nonSMG_850_oneband(0) & (data[0]['Radio_excess']==1) )  
    
    
    print "-"*dashnum
    
    maskQG = MaskAllQG_oneband(0) & Mask_nonSMG_850_oneband(0) & Mask_nonIRB(0)
    taskprint_latex( "IR-faint QG", maskQG )
    StackBin( maskQG )  

    
    print "-"*dashnum
    
    maskIRBQG = MaskAllQG_oneband(0) & Mask_nonSMG_850_oneband(0) & Mask_IRB(0)
    taskprint_latex( "IR-bright QG", maskIRBQG )
    StackBin( maskIRBQG )
    
    
    print "-"*dashnum
    print
    
    return

def Stack_test():
    
    if (SHRESHOLD==10000.0):
        print "shreshold = -    mean =",MEAN,"   sigma =",SIGMA,"/sqrt(N)"
    else:
        print "shreshold =",SHRESHOLD,"    mean =",MEAN,"   sigma =",SIGMA,"/sqrt(N)"
    
    dashnum = 85
    
    print "-"*dashnum
    print '{:<16s}{:<10s}{:<10s}{:<8s}{:^10s}{:>6s}{:>10s}{:>8s}{:>6s}'.format(
    "Groups","Mean_M*","Median_z","Number","Mean_flux","SNR","L_IR","SFR","sSFR")
    
    print "-"*dashnum
    
    #maskall = MaskAllQG(0)
    maskall = Mask_M(0) & ( data[0][Ks_error]<0.2 ) & ( data[0][Ks_error]>0 ) \
                & Mask_photoz(0) & Mask_class_star(0) & Mask_myclassQG(0)
    taskprint( "All QG", maskall )
    StackBin( maskall )
    
    print "-"*dashnum
    print
    
    return
    

def PrintStacking450():
    
    if (SHRESHOLD==10000.0):
        print "shreshold = -    mean =",MEAN,"   sigma =",SIGMA,"/sqrt(N)"
    else:
        print "shreshold =",SHRESHOLD,"    mean =",MEAN,"   sigma =",SIGMA,"/sqrt(N)"
    
    dashnum = 85
    
    print "-"*dashnum
    print '{:<16s}{:<10s}{:<10s}{:<8s}{:^10s}{:>6s}{:>10s}{:>8s}{:>6s}'.format(
    "Groups","Mean_M*","Median_z","Number","Mean_flux","SNR","L_IR","SFR","sSFR")
    
    print "-"*dashnum
    
    Mask_450Region = MaskRegion_nomask(CENTER_RA, CENTER_DEC, RADIUS)
    
    #QGs without 850um detection    #STACK_MAX,STACK_MIN = -0.1,0.1
    taskprint( "All QG",          MaskAllQG(0) & Mask_nonSMG_450(0) & Mask_450Region                                )
    
    #QGs with 24um detection        #STACK_MAX,STACK_MIN = -0.35,0.35
    taskprint( "24-micron ctp.",  MaskAllQG(0) & Mask_nonSMG_450(0) & (data[0]['24MICRON']    ==1) & Mask_450Region )
    
    #QGs with 3GHz detection        #STACK_MAX,STACK_MIN = -0.2,0.2
    taskprint( "3-GHz ctp.",      MaskAllQG(0) & Mask_nonSMG_450(0) & (data[0]['3GHZ']        ==1) & Mask_450Region )
    
    #QGs with label "SFG"           #STACK_MAX,STACK_MIN = -0.6,0.6
    taskprint( "3-GHz ctp.: SFG", MaskAllQG(0) & Mask_nonSMG_450(0) & (data[0]['Radio_excess']==0) & Mask_450Region )
    
    #QGs with label "AGN"           #STACK_MAX,STACK_MIN = -0.2,0.2
    taskprint( "3-GHz ctp.: AGN", MaskAllQG(0) & Mask_nonSMG_450(0) & (data[0]['Radio_excess']==1) & Mask_450Region )  

    maskQG = MaskAllQG(0) & Mask_nonSMG_450(0) & Mask_nonIRB(0) & Mask_450Region
    taskprint( "IR-faint QG", maskQG )
    
    maskIRBQG = MaskAllQG(0) & Mask_nonSMG_450(0) & Mask_IRB(0) & Mask_450Region
    taskprint( "IR-bright QG", maskIRBQG )
    
    
    print "-"*dashnum
    print
    
    return

def PrintStacking450OneBand():
    
    if (SHRESHOLD==10000.0):
        print "shreshold = -    mean =",MEAN,"   sigma =",SIGMA,"/sqrt(N)"
    else:
        print "shreshold =",SHRESHOLD,"    mean =",MEAN,"   sigma =",SIGMA,"/sqrt(N)"
    
    dashnum = 85
    
    print "-"*dashnum
    print '{:<16s}{:<10s}{:<10s}{:<8s}{:^10s}{:>6s}{:>10s}{:>8s}{:>6s}'.format(
    "Groups","Mean_M*","Median_z","Number","Mean_flux","SNR","L_IR","SFR","sSFR")
    
    print "-"*dashnum
    
    global LIR_iter
    LIR_iter = 31
    
    Mask_450Region = MaskRegion_nomask(CENTER_RA, CENTER_DEC, RADIUS)
    
    #QGs without 850um detection    #STACK_MAX,STACK_MIN = -0.1,0.1
    taskprint_latex( "All QG",          MaskAllQG_oneband(0) & Mask_nonSMG_450_oneband(0) & Mask_450Region                                )
    
    #QGs with 24um detection        #STACK_MAX,STACK_MIN = -0.35,0.35
    taskprint_latex( "24-micron ctp.",  MaskAllQG_oneband(0) & Mask_nonSMG_450_oneband(0) & (data[0]['24MICRON']    ==1) & Mask_450Region )
    
    #QGs with 3GHz detection        #STACK_MAX,STACK_MIN = -0.2,0.2
    taskprint_latex( "3-GHz ctp.",      MaskAllQG_oneband(0) & Mask_nonSMG_450_oneband(0) & (data[0]['3GHZ']        ==1) & Mask_450Region )
    
    #QGs with label "SFG"           #STACK_MAX,STACK_MIN = -0.6,0.6
    taskprint_latex( "3-GHz ctp.: SFG", MaskAllQG_oneband(0) & Mask_nonSMG_450_oneband(0) & (data[0]['Radio_excess']==0) & Mask_450Region )
    
    #QGs with label "AGN"           #STACK_MAX,STACK_MIN = -0.2,0.2
    taskprint_latex( "3-GHz ctp.: AGN", MaskAllQG_oneband(0) & Mask_nonSMG_450_oneband(0) & (data[0]['Radio_excess']==1) & Mask_450Region )  

    maskQG = MaskAllQG_oneband(0) & Mask_nonSMG_450_oneband(0) & Mask_nonIRB(0) & Mask_450Region
    taskprint_latex( "IR-faint QG", maskQG )
    
    maskIRBQG = MaskAllQG_oneband(0) & Mask_nonSMG_450_oneband(0) & Mask_IRB(0) & Mask_450Region
    taskprint_latex( "IR-bright QG", maskIRBQG )
    
    
    print "-"*dashnum
    print
    
    return

def PlotMan2016( axs_list ):
    
    Man_QG_y = [ [0.7,0.8,0.5,0.3],      [3.4,3.1,2.1,1.7],
                 [8.1,5.8,4.9,2.0],      [14.8,10.2,11.7,1e-10],
                 [50.1,36.3,1e-10,1e-10],[35.5,1e-10,1e-10,1e-10]    ]
    
    Man_IRBQG_y = [ [1e-10,1e-10,1e-10,1e-10], [58.9,91.2,1e-10,1e-10],
                    [89.1,85.1,89.1,1e-10],    [134.9,131.8,95.5,1e-10],
                    [269.2,218.8,1e-10,1e-10], [331.1,1e-10,1e-10,1e-10]   ]
    
    Man_SFG_y = [ [8.1,6.3,4.0,2.0],       [23.4,17.8,12.3,5.8],
                  [55.0,39.8,20.9,10.7],   [87.1,47.9,27.5,1e-10],
                  [147.9,69.2,1e-10,1e-10],[204.2,91.2,1e-10,1e-10]   ]
    
    Man_x = [11.2,10.8,10.4,10.0]
    
    for i in range(len(axs_list)):
        axs_list[i].scatter(Man_x,np.log10(Man_QG_y[i]),25,color='r',alpha=0.2,label="Man+2016: IR-faint QG")
        axs_list[i].scatter(Man_x,np.log10(Man_IRBQG_y[i]),25,color='C4',alpha=0.2,label="Man+2016: IR-bright QG",marker='D')
        axs_list[i].scatter(Man_x,np.log10(Man_SFG_y[i]),25,color='b',alpha=0.2,label="Man+2016: SFG",marker='s')
        
    return axs_list

def PlotSpeagle2014( axs_list ):
    
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    t = cosmo.age(z).value
    #t = [9.5874,6.6247,5.0223,3.6477,2.0]
    x = np.linspace(9.3, 11.7, 400)
    for i in range(len(axs_list)):
        y = (0.84 - 0.026*t[i])*x - (6.51-0.11*t[i])
        axs_list[i].plot(x,y,'k',label="Speagle+2014")
        axs_list[i].fill_between(x, y-0.9, y+0.9,color='k',alpha=0.1)
        axs_list[i].fill_between(x, y-0.6, y+0.6,color='k',alpha=0.1)
        axs_list[i].fill_between(x, y-0.3, y+0.3,color='k',alpha=0.1)
        #axs_list[i].annotate(z_range[i], xy=(0.08, 0.85), xycoords='axes fraction', fontsize=12,\
        #             xytext=(-5, 5), textcoords='offset points',\
        #             ha='left', va='bottom')
        axs_list[i].annotate(z_range[i], xy=(1, 0), xycoords='axes fraction', fontsize=12,\
                     xytext=(-5, 5), textcoords='offset points',\
                     ha='right', va='bottom')
        axs_list[i].text(9.5,(0.84 - 0.026*t[i])*9.5 - (6.51-0.11*t[i])+0.45,"z="+str(z[i]),rotation=180/np.pi*np.arctan(0.84-0.026*t[i])-20)
    
    return axs_list

def SSFR( axs_list,title, MASK, inputcolor, inputmarker ):
    
    global LIR_iter

    for i in range(len(z_list)):
        for j in range(2):
            mean_mass,SFR,SFRerr_low,SFRerr_up = taskSSFR_WH(z_list[i]&m_list[j]&MASK)
            #mean_mass,SFR,masserr_up,masserr_low = taskSSFR(z_list[i]&m_list&MASK)
                #ax_list[i].errorbar(mean_mass, SFR,xerr=np.array([[masserr_low,masserr_up]]).T,color=inputcolor,alpha=0.5)
            if(i==1)&(j==1):
                
                if (LIR_limitmask_WH[LIR_iter-1]==0):
                    axs_list[i].scatter(mean_mass,SFR,60,\
                           color=inputcolor,alpha=0.8,label=title,marker=inputmarker)
                    axs_list[i].errorbar(mean_mass, SFR,yerr=np.array([[SFRerr_low,SFRerr_up]]).T,
                           capsize=5,color=inputcolor)
            else:
                if (LIR_limitmask_WH[LIR_iter-1]==1):
                    axs_list[i].scatter(mean_mass,SFR,60,\
                           color=inputcolor,alpha=0.8,marker=inputmarker)
                    axs_list[i].errorbar(mean_mass,SFR,yerr=0.7,uplims=True,\
                           color=inputcolor,alpha=0.8)#,xerr=0.1
                    
                else:
                    axs_list[i].scatter(mean_mass,SFR,60,\
                           color=inputcolor,alpha=0.8,marker=inputmarker)
                    axs_list[i].errorbar(mean_mass, SFR,yerr=np.array([[SFRerr_low,SFRerr_up]]).T,
                           capsize=5,color=inputcolor)
    return

def PlotSSFR( MASK, FLAG ):

    ### set figure
    fig, axs = plt.subplots(2, 3,figsize=(10,6), sharex=True, sharey=True,
                        gridspec_kw={'hspace': 0, 'wspace': 0})
    if ( FLAG==0 ): axs[-1, -1].axis('off')
    
    (ax1,ax2,ax3),(ax4,ax5,ax6) = axs
    if ( FLAG==1 ): axs_list = [ax1,ax2,ax3,ax4,ax5,ax6]
    else: axs_list = [ax1,ax2,ax3,ax4,ax5]
    
    ### plot QGs
    global LIR_iter
    if ( FLAG==1 ): LIR_iter = 6
    maskQG =  MASK & Mask_nonIRB(0)
    SSFR(axs_list,'IR-radio-faint QG', maskQG,'C3','o')
    
    if ( FLAG==1 ): LIR_iter = 19
    maskIRBQG = MASK & Mask_IRB(0)
    SSFR(axs_list,'IR-radio-bright QG', maskIRBQG,'C4','D')
    
    for i in range(len(axs_list)):
        axs_list[i].set_xlim([9.3, 11.7])
        axs_list[i].set_ylim([-1.5, 3.1])
    
    '''
    #plot 450                     
                     
    SFRerr_low = np.log10( 3.4e10*1.7e-10 ) - np.log10( (3.4e10-3.4e10/4.2)*1.7e-10 )
    SFRerr_up = np.log10( (3.4e10+3.4e10/4.2)*1.7e-10 ) - np.log10( 3.4e10*1.7e-10 )
    for i in range(5):
        axs_list[i].scatter(11.0,np.log10(3.4e10*1.7e-10),60, color='C6',alpha=0.8)
        axs_list[i].errorbar(11.0, np.log10(3.4e10*1.7e-10),yerr=np.array([[SFRerr_low,SFRerr_up]]).T,
                               capsize=5,color='C6')
                        
        axs_list[i].errorbar(10.7,np.log10(6.4e8*1.7e-10),xerr=0.1,yerr=0.7,uplims=True,\
                               color='C6',alpha=0.8)
    '''
    
    ### plot Man+2016
    axs_list = PlotMan2016( axs_list )
    
    #### plot MS from Speagle+2014
    axs_list = PlotSpeagle2014( axs_list )

    
    ### plot legend
    if ( FLAG==0 ):
        axs_list[1].legend(loc='lower left', bbox_to_anchor= (1.1, -0.9), ncol=1, 
                borderaxespad=0, frameon=False)
    else: axs_list[1].legend(loc='lower left', bbox_to_anchor= (2.1, 0.2), ncol=1, 
                frameon=False,fontsize=12)

    
    for ax in axs.flat:
        ax.label_outer()       
   
    fig.text(0.5, 0.04, 'log( stellar mass($M_{\odot}$) )', ha='center', fontdict = {'fontsize' : 14})
    fig.text(0.06, 0.5, 'log( SFR($M_{\odot}$/yr) )', va='center', rotation='vertical', fontdict = {'fontsize' : 14})
    
    fig.savefig('SFRvsM.png', bbox_inches = 'tight', format='png', dpi=400)
    
    return

def PlotSSFR_fourband():
    
    ### set redshift
    global z_list, m_list, z, z_range
    z_list = [ Mask_z(0,0,0.5),Mask_z(0,0.5,1.0),Mask_z(0,1.0,1.5),Mask_z(0,1.5,2.0),Mask_z(0,2.0,100) ]
    m_list = [ Mask_mass(0,0,10.5),Mask_mass(0,10.5,100) ]
    z = [0.35,0.8,1.2,1.65,2.4]
    z_range = ["0 < z "+u"\u2264"+" 0.5","0.5 < z "+u"\u2264"+" 1.0",\
               "1.0 < z "+u"\u2264"+" 1.5","1.5 < z "+u"\u2264"+" 2.0","z > 2.0"]
    
    PlotSSFR( MaskAllQG(0) & Mask_nonSMG_850(0), 0 )
    
    return

def PlotSSFR_oneband():
    
    ### set redshift
    global z_list, m_list, z, z_range
    z_list = [ Mask_z(0,0,0.5),Mask_z(0,0.5,1.0),Mask_z(0,1.0,1.5),Mask_z(0,1.5,2.0),Mask_z(0,2.0,2.5),Mask_z(0,2.5,100) ]
    m_list = [ Mask_mass(0,0,10.5),Mask_mass(0,10.5,100) ]
    z = [0.35,0.8,1.2,1.7,2.25,2.7]
    z_range = ["0 < $z$ "  +u"\u2264"+" 0.5",  "0.5 < $z$ "+u"\u2264"+" 1.0",\
               "1.0 < $z$ "+u"\u2264"+" 1.5",  "1.5 < $z$ "+u"\u2264"+" 2.0",\
               "2.0 < $z$ "+u"\u2264"+" 2.5",  "$z$ > 2.5"]
    
    PlotSSFR( MaskAllQG_oneband(0) & Mask_nonSMG_850_oneband(0), 1 )
    
    return
    
def PlotSSFR850450():
    plt.figure(figsize=(8,6))
        
    #plot MS from Speagle+2014
    z = 0.9
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    t = cosmo.age(z).value
    x = np.linspace(9.8, 11.7, 400)
    y = (0.84 - 0.026*t)*x - (6.51-0.11*t)
    plt.plot(x,y,'k',label="Speagle+2014")
    plt.fill_between(x, y-0.9, y+0.9,color='k',alpha=0.1)
    plt.fill_between(x, y-0.6, y+0.6,color='k',alpha=0.1)
    plt.fill_between(x, y-0.3, y+0.3,color='k',alpha=0.1)
    plt.text(10.6,(0.84 - 0.026*t)*10.6 - (6.51-0.11*t+0.1),"z=0.9",
             rotation=180/np.pi*np.arctan(0.84-0.026*t)-20)
    
    label_list = ["All QG candidates","24um ctp.","3GHz ctp.",
                  "3GHz ctp.:SFG","3GHz ctp.:AGN","IR-bright QG","QG"]
    marker_list = ['o','s','D','d','x','^','v']
    
    mass_list_850 = [10.7, 11.0, 11.1, 11.1, 11.2, 11.0, 10.7]
    LIR_list_850 = np.array([1.0e10,1.1e11,6.8e10,1.4e11,1.7e9,1.2e11,6.7e9])
    LIR_limitmask_850 = [0,0,0,0,1,0,0]
    LIR_errlen_850 = np.array([1.0e10/6.2,1.1e11/8.0,6.8e10/3.6,1.4e11/5.1,
                                8.5e9,
                                1.2e11/8.4,6.7e9/3.7])
    
    mass_list_450 = [10.7,10.9,11.2,11.1,11.2,11.0,10.7]
    LIR_list_450 = np.array([2.6e9,3.4e10,3.5e10,6.5e10,1.5e10,3.4e10,6.4e8])
    LIR_limitmask_450 = [1,0,0,0,1,0,1]
    LIR_errlen_450 = np.array([1.6e9,
                               3.4e10/3.9,3.5e10/3.0,6.5e10/2.8,
                               9.7e9,
                               3.4e10/4.2,
                               1.6e9])
    
    #legend management
    plt.scatter(-1000,-1000,60,color='C6',alpha=0.8,label='450um,SNR'+u"\u2265"+'2.0')
    plt.scatter(-1000,-1000,60,color='C4',alpha=0.8,label='850um,SNR'+u"\u2265"+'2.0')
    plt.scatter(-1000,-1000,60,color='C6',alpha=0.2,label='450um,SNR<2.0')
    plt.scatter(-1000,-1000,60,color='C4',alpha=0.2,label='850um,SNR<2.0')
    for i in range( len(label_list) ):
        plt.scatter(-1000,-1000,color='k',alpha=0.2,marker=marker_list[i],label=label_list[i])
    
    #plot the points
    for i in range( len(label_list) ):
        
        SFR_450 = np.log10(LIR_list_450[i]*1.7e-10)
        SFRerr_up_450 = np.log10((LIR_list_450[i]+LIR_errlen_450[i])*1.7e-10) - SFR_450
        if (LIR_list_450[i]-LIR_errlen_450[i])<=0:
            SFRerr_low_450 = 100
        else:
            SFRerr_low_450 = SFR_450 - np.log10((LIR_list_450[i]-LIR_errlen_450[i])*1.7e-10)
            
        inputalpha = 0.0
        if LIR_limitmask_450[i]==0:
            inputalpha = 0.8
        else:
            inputalpha = 0.2
        
        plt.scatter(mass_list_450[i],SFR_450,60,
                    marker=marker_list[i],color='C6',alpha=inputalpha)
        plt.errorbar(mass_list_450[i], SFR_450,
                     yerr=np.array([[SFRerr_low_450,SFRerr_up_450]]).T,
                     capsize=5,color='C6',alpha=inputalpha)
 
        SFR_850 = np.log10(LIR_list_850[i]*1.7e-10)
        SFRerr_up_850 = np.log10((LIR_list_850[i]+LIR_errlen_850[i])*1.7e-10) - SFR_850
        if (LIR_list_850[i]-LIR_errlen_850[i])<=0:
            SFRerr_low_850 = 100
        else:
            SFRerr_low_850 = SFR_850 - np.log10((LIR_list_850[i]-LIR_errlen_850[i])*1.7e-10)
            
        inputalpha = 0.0
        if LIR_limitmask_850[i]==0:
            inputalpha = 0.8
        else:
            inputalpha = 0.2
        
        plt.scatter(mass_list_850[i],SFR_850,60,
                    marker=marker_list[i],color='C4',alpha=inputalpha)
        plt.errorbar(mass_list_850[i], SFR_850,
                     yerr=np.array([[SFRerr_low_850,SFRerr_up_850]]).T,
                     capsize=5,color='C4',alpha=inputalpha)

    
    
    plt.axis([10.2,11.4,-1.3,1.8])
    plt.legend()
    plt.xlabel('log( stellar mass($M_{\odot}$) )', fontdict = {'fontsize' : 14})
    plt.ylabel('log( SFR($M_{\odot}$/yr) )', fontdict = {'fontsize' : 14})
    
    return
    
def PlotSSFR850450_tick():
    
    LIR_list_850 = np.array([1.0e10,1.1e11,6.8e10,1.4e11,1.7e9,1.2e11,6.7e9])
    LIR_limitmask_850 = [0,0,0,0,1,0,0]
    LIR_errlen_850 = np.array([1.0e10/6.2,1.1e11/8.0,6.8e10/3.6,1.4e11/5.1,
                                8.5e9,
                                1.2e11/8.4,6.7e9/3.7])
    
    LIR_list_450 = np.array([2.6e9,3.4e10,3.5e10,6.5e10,1.5e10,3.4e10,6.4e8])
    LIR_limitmask_450 = [1,0,0,0,1,0,1]
    LIR_errlen_450 = np.array([1.6e9,
                               3.4e10/3.9,3.5e10/3.0,6.5e10/2.8,
                               9.7e9,
                               3.4e10/4.2,
                               1.6e9])
        
    fig, axes = plt.subplots(figsize=(10,6)) 
    plt.setp(axes, xticks=[y_axes+1 for y_axes in range( len(LIR_list_850) )],
             xticklabels=["All QG","24um ctp.","3GHz ctp.",
                          "3GHz ctp.:\nSFG","3GHz ctp.:\nAGN","IR-bright\nQG","QG"])
        
    for i in range( len(LIR_list_850) ):
        
        SFR_450 = np.log10(LIR_list_450[i]*1.7e-10)
        SFRerr_up_450 = np.log10((LIR_list_450[i]+LIR_errlen_450[i])*1.7e-10) - SFR_450
        if (LIR_list_450[i]-LIR_errlen_450[i])<=0:
            SFRerr_low_450 = 100
        else:
            SFRerr_low_450 = SFR_450 - np.log10((LIR_list_450[i]-LIR_errlen_450[i])*1.7e-10)
            
        inputalpha = 0.0
        if LIR_limitmask_450[i]==0: inputalpha = 0.8
        else: inputalpha = 0.2
            
        plt.errorbar(i+1,SFR_450,yerr=np.array([[SFRerr_low_450,SFRerr_up_450]]).T,
                     color='C6', fmt='o', capsize=5, alpha=inputalpha)
        
        SFR_850 = np.log10(LIR_list_850[i]*1.7e-10)
        SFRerr_up_850 = np.log10((LIR_list_850[i]+LIR_errlen_850[i])*1.7e-10) - SFR_850
        if (LIR_list_850[i]-LIR_errlen_850[i])<=0:
            SFRerr_low_850 = 100
        else:
            SFRerr_low_850 = SFR_850 - np.log10((LIR_list_850[i]-LIR_errlen_850[i])*1.7e-10)
            
        inputalpha = 0.0
        if LIR_limitmask_850[i]==0: inputalpha = 0.8
        else: inputalpha = 0.2
            
        plt.errorbar(i+1,SFR_850,yerr=np.array([[SFRerr_low_850,SFRerr_up_850]]).T,
                     color='C4', fmt='o', capsize=5, alpha=inputalpha)

    plt.axis([0,8,-1.3,1.8])
    
    #legend management
    plt.scatter(-1000,-1000,60,color='C6',alpha=0.8,label='450um,SNR'+u"\u2265"+'2.0')
    plt.scatter(-1000,-1000,60,color='C4',alpha=0.8,label='850um,SNR'+u"\u2265"+'2.0')
    plt.scatter(-1000,-1000,60,color='C6',alpha=0.2,label='450um,SNR<2.0')
    plt.scatter(-1000,-1000,60,color='C4',alpha=0.2,label='850um,SNR<2.0')
    plt.legend()
    
    plt.ylabel('log( SFR($M_{\odot}$/yr) )', fontdict = {'fontsize' : 14})
    plt.show()
    
    return
    
    
def PlotSSFR850450_tick2():
    
    LIR_list_850 = np.array([1.0e10,1.1e11,6.8e10,1.4e11,1.7e9,1.2e11,6.7e9])
    LIR_limitmask_850 = [0,0,0,0,1,0,0]
    LIR_errlen_850 = np.array([1.0e10/6.2,1.1e11/8.0,6.8e10/3.6,1.4e11/5.1,
                                8.5e9,
                                1.2e11/8.4,6.7e9/3.7])
    
    LIR_list_450 = np.array([2.6e9,3.4e10,3.5e10,6.5e10,1.5e10,3.4e10,6.4e8])
    LIR_limitmask_450 = [1,0,0,0,1,0,1]
    LIR_errlen_450 = np.array([1.6e9,
                               3.4e10/3.9,3.5e10/3.0,6.5e10/2.8,
                               9.7e9,
                               3.4e10/4.2,
                               1.6e9])
        
    fig, axes = plt.subplots(figsize=(10,6)) 
    plt.setp(axes, xticks=[y_axes+1 for y_axes in range( len(LIR_list_850) )],
             xticklabels=["All QG","24um ctp.","3GHz ctp.",
                          "3GHz ctp.:\nSFG","3GHz ctp.:\nAGN","IR-bright\nQG","QG"])
        
    for i in range( len(LIR_list_850) ):
        
        inputalpha = 0.0
        if LIR_limitmask_450[i]==0: inputalpha = 0.8
        else: inputalpha = 0.2
            
        plt.errorbar(i+1,LIR_list_450[i]*1.7e-10,yerr=LIR_errlen_450[i]*1.7e-10,
                     color='C6', fmt='o', capsize=5, alpha=inputalpha)
        
            
        inputalpha = 0.0
        if LIR_limitmask_850[i]==0: inputalpha = 0.8
        else: inputalpha = 0.2
            
        plt.errorbar(i+1,LIR_list_850[i]*1.7e-10,yerr=LIR_errlen_850[i]*1.7e-10,
                     color='C4', fmt='o', capsize=5, alpha=inputalpha)

    plt.axis([0,8,0,30])
    
    #legend management
    plt.scatter(-1000,-1000,60,color='C6',alpha=0.8,label='450um,SNR'+u"\u2265"+'2.0')
    plt.scatter(-1000,-1000,60,color='C4',alpha=0.8,label='850um,SNR'+u"\u2265"+'2.0')
    plt.scatter(-1000,-1000,60,color='C6',alpha=0.2,label='450um,SNR<2.0')
    plt.scatter(-1000,-1000,60,color='C4',alpha=0.2,label='850um,SNR<2.0')
    plt.legend()
    
    plt.ylabel('SFR($M_{\odot}$/yr)', fontdict = {'fontsize' : 14})
    plt.show()
    
    return

def PlotLensed():
    
    ###read file
    hdu_list = fits.open("/Users/yuhsuan/Desktop/COSMOS.Ks.original_psf.v5.fits")
    #hdu_list.info()
    image_data = hdu_list[0].data
    wcs = WCS(hdu_list[0].header).celestial
    
    ###chop the image
    x_center, y_center = wcs.all_world2pix([[150.1001356, 2.2971499]], 0)[0]
    #print x_center, y_center
    y_min, y_max = int(y_center)-15, int(y_center)+16
    x_min, x_max = int(x_center)-15, int(x_center)+16

    imagebox = np.copy( image_data[y_min:y_max:1,x_min:x_max:1] )
    
    
    fig, ax = plt.subplots()
    plt.imshow(imagebox, cmap='gray', vmin=-2200, vmax=2500)
    #plt.colorbar()
    plt.gca().invert_yaxis()
    
    ###plot QG sample
    circle1 = plt.Circle( (x_center-x_min,y_center-y_min), 6.6666869, edgecolor="k", facecolor='none' )
    ax.add_artist(circle1)
    
    ###plot ALMA sources
    a1_x, a1_y = wcs.all_world2pix([[150.0998696, 2.2972549]], 0)[0]
    a2_x, a2_y = wcs.all_world2pix([[150.1001325, 2.2970430]], 0)[0]
    circle_a1 = plt.Circle( (a1_x-x_min,a1_y-y_min), 3.3333435, edgecolor="pink", facecolor='none' )
    circle_a2 = plt.Circle( (a2_x-x_min,a2_y-y_min), 3.3333435, edgecolor="pink", facecolor='none' )
    ax.add_artist(circle_a1)
    ax.add_artist(circle_a2)
    plt.text( a1_x-x_min-2.5, a1_y-y_min-5, "AS2COS0005.1", fontsize=12, color="pink" )
    plt.text( a2_x-x_min-2.5, a2_y-y_min-5, "AS2COS0005.2", fontsize=12, color="pink" )
    
    ###plotting scale
    plt.plot( (1.0,1.0+6.6666869), (1.0,1.0), 'k')
    plt.text( 1.0+6.6666869*0.07, 1.4, "1 arcsec", fontsize=12, color="k" )
    
    fig.savefig('lensed_Ks.png', bbox_inches = 'tight', format='png', dpi=400)
    
    return

# =============================================================================
# main code
# =============================================================================

time1 = time.time()

###set columns in the main catalog
color1     = "MNUV"
color2     = "MR"
color3     = "MJ"
color      = [ color1, color2, color3 ]
photoz     = "REDSHIFT"
V_error    = "V_MAGERR_APER3"
ip_error   = "ip_MAGERR_APER3"
J_error    = "J_MAGERR_APER3"
Ks_error   = "Ks_MAGERR_APER3"
Ks_mag     = "Ks_MAG_APER3"
mass       = "MASS_MED"
class_star = "TYPE"
class_SFG  = "CLASS"
ra         = "ALPHA_J2000"
dec        = "DELTA_J2000"

mass_error_upper = "MASS_MED_MAX68"
mass_error_lower = "MASS_MED_MIN68"


ra         = [None]*1
dec        = [None]*1
ra[0]      = "ALPHA_J2000"
dec[0]     = "DELTA_J2000"

###read catalog
catalog    = [None]*1
catalog[0] = 'COSMOS2015_merged.fits' 
data       = [None]*1
x,y        = [],[]
ReadCatalog( 0,catalog[0] )

###read image file
path  = "/Users/yuhsuan/Documents/research/05WH/data/COSMOS/"
IMAGE = path+'04_COSMOS450_850/S2COSMOS/maps/S2COSMOS_20180927_850_fcf_mf_crop.fits'
RMS   = path+'04_COSMOS450_850/S2COSMOS/maps/S2COSMOS_20180927_850_err_mf_crop.fits'
image_data, wcs = ReadImage(IMAGE,1)
rms_data        = ReadImage(RMS,0)

###set shreshold
SHRESHOLD = 3.0

######### random stacking
'''
print 'random stacking'

MASK2 = Mask_M(0) & Mask_photoz(0) & Mask_error(1,0.1,0) & Mask_class_star(0)
MASK3 = Mask_myclassQG(0)
MASK = MASK2 & MASK3

RandomStackingTask()
'''

######### stack QG

if   ( SHRESHOLD==10000.0 ):  MEAN, SIGMA =  0.000, 1.326
elif ( SHRESHOLD==3.5     ):  MEAN, SIGMA = -0.029, 1.265
elif ( SHRESHOLD==3.0     ):  MEAN, SIGMA = -0.044, 1.279
elif ( SHRESHOLD==2.0     ):  MEAN, SIGMA = -0.130, 1.166


IMAGE_MAX, IMAGE_MIN = -10,10


sigma_plot = 1.0

#four-band
'''
LIR_WH = [1.4e9*sigma_plot, 6.8e9*sigma_plot/1.3,
          1.6e9*sigma_plot/0.5,8.0e9,
          7.6e9*sigma_plot,3.6e10,
          3.5e10*sigma_plot/1.3,3.7e10,
          6.9e10*sigma_plot/0.6,1.2e11,
          
          6.6e10,2.0e10*sigma_plot/1.6,
          1.0e11,3.2e10,
          1.7e11*sigma_plot/1.6,2.4e11,
          1.1e12,1.0e12,
          1.1e12*sigma_plot/1.6,3.9e11]

#mean_flux scale to 1 sigma + sigma_plot
#except for zero or minus mean_flux change to 0+1 sigma
#LIR_errlen_WH = [1.4e9,6.8e9/1.3,
#                 1.6e9/0.5,8.0e9/2.2,
#                 7.6e9,3.6e10/3.5,
#                 3.5e10/1.3,3.7e10/2.1,
#                 6.9e10/0.6,1.2e11/3.4,
#                 
#                 6.6e10/2.0,2.0e10/1.6,
#                 1.0e11/2.8,3.2e10/2.6,
#                 1.7e11/1.6,2.4e11/5.2,
#                 1.1e12/3.0,1.0e12/6.4,
#                 1.1e12/1.6,3.9e11/2.8]

#mean_flux scale to 1 sigma + sigma_plot
#except for SNR<2.0 change to 0+1 sigma
LIR_errlen_WH = [1.4e9,4.9e9,
                 3.2e9,8.0e9/2.2,
                 7.6e9,3.6e10/3.5,
                 3.2e10,3.7e10/2.1,
                 1.1e11,1.2e11/3.4,
                 
                 6.6e10/2.0,9.0e9,
                 1.0e11/2.8,3.2e10/2.6,
                 9.7e10,2.4e11/5.2,
                 1.1e12/3.0,1.0e12/6.4,
                 6.5e11,3.9e11/2.8]

for i in range( len(LIR_errlen_WH) ):
    LIR_errlen_WH[i] = LIR_errlen_WH[i] * sigma_plot

LIR_limitmask_WH = [1,1,1,0,1,0,1,0,1,0,   0,1,0,0,1,0,0,0,1,0]
LIR_iter = 0
'''


LIR_iter = 0

SNR = [ 6.3, 8.2, 3.5, 5.0, 0.1,
        3.8,
        -0.4,1.4, 0.4, 2.2,
        -0.5,3.4, 1.3, 2.3,
        0.9, 3.2, -0.3,1.1,
        8.6,
        2.0, 1.6, 2.8, 2.7,
        2.0, 5.1, 2.5, 6.6,
        2.4, 2.4, 0.5, 1.8 ]

LIR_limitmask_WH = []

for i in range( len(SNR) ):
    if ( SNR[i] < 2.0 ): LIR_limitmask_WH.append(1)
    else: LIR_limitmask_WH.append(0)

LIR_WH = [ 1.0393994e+10, 1.1675221e+11, 6.7596636e+10, 1.3854581e+11, 1.6799213e+09,
           6.7196851e+09,
           2.1875100e+09, 6.8162994e+09, 1.5864483e+09, 8.0078843e+09,
           5.6807255e+09, 3.5750400e+10, 2.0815037e+10, 3.7364623e+10,
           2.9396089e+10, 1.2770901e+11, 9.9575061e+10, 6.0067898e+10,
           1.1661365e+11,
           6.8301539e+10, 2.0047333e+10, 1.0498772e+11, 3.0109917e+10,
           2.1107203e+11, 2.2635037e+11, 8.0875008e+11, 1.0234224e+12,
           1.2881062e+12, 3.8623161e+11, 2.8080578e+11, 3.0794453e+11,
           1.0393994e+10, 8.0367031e+11, 6.6293680e+11, 1.0435518e+12, 1.9249000e+11,
           8.5229484e+09, 7.6447669e+11 ]
'''
for i in range( len(SNR) ):
    if ( SNR[i] < 2.0 )&( SNR[i] > 0.0 ):
        LIR_WH[i] *= sigma_plot
        LIR_WH[i] /= SNR[i]
'''

LIR_errlen_WH = []

for i in range( len(SNR) ):
    if ( SNR[i] >= 2.0 ): LIR_errlen_WH.append( LIR_WH[i]*sigma_plot/SNR[i] )
    else: LIR_errlen_WH.append( -99.9 )




#PrintStacking()
#PrintStackingOneBand()
#Stack_test()

#PlotSSFR_fourband()
PlotSSFR_oneband()

#PlotSSFR850450()
#PlotSSFR850450_tick()
#PlotSSFR850450_tick2()


#PlotLensed()


###### 450

###read image file
path  = "/Users/yuhsuan/Documents/research/05WH/data/COSMOS/"
IMAGE = path+'04_COSMOS450_850/STUDIES/STUDIES+S2CLS+Casey_450_mf.fits'
RMS   = path+'04_COSMOS450_850/STUDIES/STUDIES+S2CLS+Casey_450_mf_rms.fits'
image_data, wcs = ReadImage(IMAGE,1)
rms_data        = ReadImage(RMS,0)

###set shreshold
SHRESHOLD = 3.0

if   ( SHRESHOLD==10000.0 ):  MEAN, SIGMA =  0.000, 2.414
elif ( SHRESHOLD==3.5     ):  MEAN, SIGMA = -0.126, 2.073
elif ( SHRESHOLD==3.0     ):  MEAN, SIGMA = -0.153, 2.005
elif ( SHRESHOLD==2.0     ):  MEAN, SIGMA = -0.290, 1.931

CENTER_RA,CENTER_DEC = '10h00m25.0s', '2d24m22.0s'
RADIUS               = 0.2*u.degree

#PrintStacking450()
#PrintStacking450OneBand()

time2 = time.time()
print 'done! time :', time2-time1 , 'sec'