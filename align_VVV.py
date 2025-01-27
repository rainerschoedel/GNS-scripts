#!/usr/bin/env python
# coding: utf-8

# In[13]:


import numpy as np
import matplotlib.pyplot as plt
import astroalign as aa
from astropy.io import fits
from scipy.spatial import distance
from skimage.transform import warp
from astropy.table import Table
import pandas as pd
import sys

# VVV --> HAWKI
scale = 0.34/0.106

# number of brightest stars in te lists that wil be used
# (this will speed up things a lot)
n_bright_hawki = 100
n_bright_vvv = 200


#get parameters from script call
field = str(sys.argv[1])
band  = sys.argv[2]
chip  = str(sys.argv[3])


#define paths
VVV_path  ='/home/data/VVV/Fields/J/'
lnx_path  = f'/home/data/GNS/2015/{band}/{field}/ims/'
lnx_stars = f'/home/data/GNS/2015/{band}/{field}/data/'
data_path = f'/home/data/GNS/2015/{band}/{field}/data/'
cube_path = f'/home/data/GNS/2015/{band}/{field}/cubes/'

#Read VVV data
img_vvv = fits.getdata(VVV_path + f'Field{field}.fits.gz')
vvv_img = pd.DataFrame(np.array(img_vvv).byteswap().newbyteorder())
xsize = vvv_img.shape[1]
ysize = vvv_img.shape[0]
xh = round(xsize/2)
yh = round(ysize/2)

print(f'Size of VVV image is {xsize} x {ysize}.')
tab_vvv = Table.read(VVV_path + f'Field{field}_stars.txt',format='ascii')
n_vvv = len(tab_vvv)
x_vvv = tab_vvv['x']
y_vvv = tab_vvv['y']
m_vvv = tab_vvv['Ks']
print(f'There are {n_vvv} stars in the VVV image.') 

#crop list for each chip
if (chip == '1'):
    idx = np.nonzero((x_vvv < xh) & (y_vvv < yh))
elif (chip == '2'):
    idx = np.nonzero((x_vvv > xh) & (y_vvv < yh))
elif (chip == '3'):
    idx = np.nonzero((x_vvv > xh) & (y_vvv > yh))
elif (chip == '4'):
    idx = np.nonzero((x_vvv < xh) & (y_vvv > yh))
x_vvv = x_vvv[idx]
y_vvv = y_vvv[idx]
print(f'There are {len(x_vvv)} stars in the corresponding section of the VVV image.')
print(f'Using the {n_bright_hawki} brightest stars in HAWK-I for alignment.')
print(f'Using the {n_bright_vvv} brightest stars in VVV for alignment.')
print(f'Brightest used star from VVV has J = {m_vvv[0]} and faintest one {m_vvv[n_bright_vvv]}.')
x_vvv = x_vvv[0:n_bright_vvv]
y_vvv = y_vvv[0:n_bright_vvv]
vvv_xy = np.array([[x, y] for x, y in zip(x_vvv, y_vvv)], dtype="float64")


#Read HAWK-I data
#limit number of stars to the brightest n_brgiht
#The input list is supposed to be ordered by descending brightness
tab_img = Table.read(data_path + 'stars_' + chip + '.txt',format='ascii')
x_img = tab_img['x'][0:n_bright_hawki]
y_img = tab_img['y'][0:n_bright_hawki]
img_xy = np.array([[x, y] for x, y in zip(x_img, y_img)], dtype="float64")

#TRANSFORMAITON INTO VVV
#save ordered lists of common stars
#===================================
#Determine transformation into VVV orig
p, (pos_img, pos_img_t) = aa.find_transform(img_xy, vvv_xy, max_control_points=200)
print("\nTranformation matrix:\n{}".format(p.params))
np.savetxt(data_path + f'aa_stars_hawki_'+chip+'.txt',pos_img)
np.savetxt(data_path + f'aa_stars_vvv_'+chip+'.txt',pos_img_t)
