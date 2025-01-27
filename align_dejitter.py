#!/usr/bin/env python
# coding: utf-8

# In[13]:


import numpy as np
import astroalign as aa
from astropy.io import fits
import sys
from astropy.table import Table

#get parameters from script call
field = sys.argv[1]
band = sys.argv[2]
n_exp = int(sys.argv[3])
chip = '1'

ref_cube = '1' #alignment will be done with respect to this exposure


cube_path = '/home/data/GNS/2015/' + band + '/' + field + '/ims/'
data_path =  '/home/data/GNS/2015/'+ band + '/' + field + '/data/'
out_path =  '/home/data/GNS/2015/' + band + '/' + field + '/cubes/'
tmp_path =  '/home/data/GNS/2015/' + band + '/' + field + '/tmp/'


cube_ref = fits.getdata(cube_path + 'chip' + chip + '_cube' + ref_cube + '.fits.gz')

tab_ref = Table.read(data_path + 'dejitter_stars_' + chip + '_' + ref_cube + '.txt',format='ascii')
x_ref = tab_ref['x']
y_ref = tab_ref['y']
ref_xy = np.array([[x, y] for x, y in zip(x_ref, y_ref)], dtype="float64")



# Now determine transformations for each offset and align the frames


for e in range(n_exp):

    tab = Table.read(data_path + 'dejitter_stars_' + chip + '_' + str(e+1) + '.txt',format='ascii')
    x_img = tab['x']
    y_img = tab['y']
    img_xy = np.array([[x, y] for x, y in zip(x_img, y_img)], dtype="float64")
    p, (pos_img, pos_img_t) = aa.find_transform(img_xy, ref_xy,max_control_points=200)
    print("\nTranformation matrix:\n{}".format(p.params))
    np.savetxt(data_path + 'dejitter_common_' + chip+'_' + str(e+1) + '.txt',pos_img,fmt='%13.3f')
    np.savetxt(data_path + 'dejitter_common_ref_' + chip +'_' + str(e+1) + '.txt',pos_img_t,fmt='%13.3f')


    print('Aligned frames of chip ' +chip +  ', exposure ' + str(e+1) + '.')

