PRO fullbpm, common_path, tmp_path, bpm_name, sky_name, dark_name, mask_name, flat_name, SIGMA_DEV_SKY = sigma_dev_sky, SIGMA_DEV_DARK = sigma_dev_dark

nax1 = 4096
nax2 = 1536

band = 'Ks'
field = '20'

common_path = '/home/data/GNS/2015/' + band + '/' + field +'/ims/'
tmp_path = '/home/data/GNS/2015/' + band + '/' + field +'/tmp/'
bpm_name = 'bpm_' + band + '.fits'
sky_name = 'sky.fits'
dark_name = 'dark.fits'
mask_name = 'mask_' + band + '.fits'
flat_name = 'flat_' + band + '.fits'
sigma_dev_sky = 5. ; sigma threshold for determining bad pixels in sky
sigma_dev_dark = 20. ; sigma threshold for determining bad pixels in dark

if not(KEYWORD_SET(sigma_dev_sky)) then sigma_dev_sky = 5.
if not(KEYWORD_SET(sigma_dev_dark)) then sigma_dev_dark = 20.

bpm = readfits(common_path + bpm_name)
mask = readfits(common_path + mask_name)
sky = readfits(common_path + sky_name)
dark = readfits(common_path + dark_name)
flat = readfits(common_path + flat_name)

; initialize dark and sky bpms
skybpm = sky
skybpm[*,*] = 0
darkbpm = dark
darkbpm[*,*] = 0

; get rid of flat  variations in sky
good = where(mask gt 0 and bpm lt 1)
sky[good] = sky[good]/flat[good]
;sky = sky - dark

; 1) Create bad pixel map from mean and standard deviation of sky in each
; quadrant
; ----------------------------------------------------------------------
  q = sky[0:2047,0:767]
 RESISTANT_Mean,q,3.0,Mean,Sigma,Num_Rej, goodvec=goodvec
 Sigma = Sigma * sqrt(n_elements(goodvec)) ; we are interested in sigma, not error of the mean
 bad = where(abs(q - Mean) gt sigma_dev_sky * Sigma, n_bad)
 q[*,*] = 0
 q[bad] = 1
 skybpm[0:2047,0:767] = q

 q = sky[2048:4095,0:767]
 RESISTANT_Mean,q,3.0,Mean,Sigma,Num_Rej, goodvec=goodvec
 Sigma = Sigma * sqrt(n_elements(goodvec)) ; we are interested in sigma, not error of the mean
 bad = where(abs(q - Mean) gt sigma_dev_sky * Sigma, n_bad)
 q[*,*] = 0
 q[bad] = 1
 skybpm[2048:4095,0:767] = q

 q = sky[2048:4095,768:1535]
 RESISTANT_Mean,q,3.0,Mean,Sigma,Num_Rej, goodvec=goodvec
 Sigma = Sigma * sqrt(n_elements(goodvec)) ; we are interested in sigma, not error of the mean
 bad = where(abs(q - Mean) gt sigma_dev_sky * Sigma, n_bad)
 q[*,*] = 0
 q[bad] = 1
 skybpm[2048:4095,768:1535] = q

 q = sky[0:2047,768:1535]
 RESISTANT_Mean,q,3.0,Mean,Sigma,Num_Rej, goodvec=goodvec
 Sigma = Sigma * sqrt(n_elements(goodvec)) ; we are interested in sigma, not error of the mean
 bad = where(abs(q - Mean) gt sigma_dev_sky * Sigma, n_bad)
 q[*,*] = 0
 q[bad] = 1
 skybpm[0:2047,768:1535] = q

 writefits, tmp_path + 'sky_bpm.fits', skybpm

; 2) Create bad pixel map from mean and standard deviation of dark in each
; quadrant
; ----------------------------------------------------------------------

 
 q = dark[0:2047,0:767]
 RESISTANT_Mean,q,3.0,Mean,Sigma,Num_Rej, goodvec=goodvec
 Sigma = Sigma * sqrt(n_elements(goodvec)) ; we are interested in sigma, not error of the mean
 bad = where(abs(q - Mean) gt sigma_dev_dark * Sigma, n_bad)
 q[*,*] = 0
 q[bad] = 1
 darkbpm[0:2047,0:767] = q

 q = dark[2048:4095,0:767]
 RESISTANT_Mean,q,3.0,Mean,Sigma,Num_Rej, goodvec=goodvec
 Sigma = Sigma * sqrt(n_elements(goodvec)) ; we are interested in sigma, not error of the mean
 bad = where(abs(q - Mean) gt sigma_dev_dark * Sigma, n_bad)
 q[*,*] = 0
 q[bad] = 1
 darkbpm[2048:4095,0:767] = q

 q = dark[2048:4095,768:1535]
 RESISTANT_Mean,q,3.0,Mean,Sigma,Num_Rej, goodvec=goodvec
 Sigma = Sigma * sqrt(n_elements(goodvec)) ; we are interested in sigma, not error of the mean
 bad = where(abs(q - Mean) gt sigma_dev_dark * Sigma, n_bad)
 q[*,*] = 0
 q[bad] = 1
 darkbpm[2048:4095,768:1535] = q

 q = dark[0:2047,768:1535]
 RESISTANT_Mean,q,3.0,Mean,Sigma,Num_Rej, goodvec=goodvec
 Sigma = Sigma * sqrt(n_elements(goodvec)) ; we are interested in sigma, not error of the mean
 bad = where(abs(q - Mean) gt sigma_dev_dark * Sigma, n_bad)
 q[*,*] = 0
 q[bad] = 1
 darkbpm[0:2047,768:1535] = q
 writefits, tmp_path + 'dark_bpm.fits', darkbpm



; Create a full dead pixel map
new = where(skybpm gt 0)
bpm[new] = 1
new = where(darkbpm gt 0)
bpm[new] = 1
; apply mask to final bpm
; to avoid huge interpolation overheads in data reduction
writefits, common_path + 'bpm.fits', bpm * mask
writefits, common_path + 'dark_bpm.fits', darkbpm * mask

print, 'fullbpm.pro ended'

END
