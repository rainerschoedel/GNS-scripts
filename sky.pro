PRO sky, field, band, NDIT

; There are only 5 jitter positions for the sky and there is serious
; contamination by stars
; All pixels more than sigma_dev * sigma above the mean of a quadrant
; will be masked as stars. These pixels will then be filled with
; random values around the Mean.
; Finally, a median sky is calculated.
;
; The mean value of the sky can vary between exposures. Therefore,
; each exposure if normalized by its mean.
; 
; The result is a normalised sky, that has to be fitted to 
; each given science exposure by determining the sky offset 
; in each science exposure.

NDIT = 20 ; number of sub-integrations
band = 'Ks'
field = '20'
all = 0   ; average individual burst exposures (all = 0) or use mean image at the end of cube (all = 1)

in_path = '/home/data/raw/2015/' + band +'/Sky/' + field + '/'
common_path = '/home/data/GNS/2015/' + band + '/' + field +'/ims/'
out_path = '/home/data/GNS/2015/' + band + '/' + field +'/ims/'
tmp_path = '/home/data/GNS/2015/' + band + '/' + field +'/tmp/'
flat_name = 'flat_' + band + '.fits'
bpm_name = 'bpm_' + band + '.fits'
mask_name = 'mask_' + band + '.fits'
dark_name = 'dark.fits'
sigma_dev = 5. ; sigma threshold for determining valid pixels


if not(KEYWORD_SET(sigma_dev)) then sigma_dev = 10.

nax1 = 4096
nax2 = 1536

dark = readfits(common_path + dark_name)
bpm = readfits(common_path + bpm_name)
mask = readfits(common_path + mask_name)
flat = readfits(common_path + flat_name)

; identify bad pixels for interpolation
bad = where(bpm gt 0 and mask gt 0)
xy = array_indices(bpm,bad)
xbad = xy[0,*]
ybad = xy[1,*]

q1_mask = mask[0:2047,0:767]
q2_mask = mask[2048:4095,0:767]
q3_mask = mask[2048:4095,768:1535]
q4_mask = mask[0:2047,768:1535]
q1_dark = dark[0:2047,0:767]
q2_dark = dark[2048:4095,0:767]
q3_dark = dark[2048:4095,768:1535]
q4_dark = dark[0:2047,768:1535]
q1_bpm = bpm[0:2047,0:767]
q2_bpm = bpm[2048:4095,0:767]
q3_bpm = bpm[2048:4095,768:1535]
q4_bpm = bpm[0:2047,768:1535]

name = ''
inlist = in_path+'list.txt'
readcol, inlist, list, FORMAT='(A)'
ncubes = n_elements(list)
print, list

skycube = fltarr(nax1,nax2,ncubes)

for ic = 0L, ncubes-1 do begin  ; start loop over frames in this cube
 nam = list[ic]
 
 ; Assume that mean of cube is contained in the last slice.
 ; This is standard in HAWK-I FASTPHOT mode.
 
 
 if all eq 1 then begin
  im = readfits(in_path+nam, NSLICE = NDIT) 
 endif else begin
  ; Compute the mean from individual images
  ; better option because first frame is often of worse quality
  ; exclude first frame and also the last one, which is the mean
  ; of all othr frames
  im_tmp = readfits(in_path+nam) 
  im_tmp = im_tmp[*,*,1:NDIT-2]
  im = mean(im_tmp,DIMENSION=3)
 endelse
;writefits, tmp_path + 'skycube.fits', im_tmp


 ; compute mean skies for each quadrant
 ; fill in bad pixels or pixels affected by stars with a random value
 ; save gain values for each quadrant
 gains = fltarr(4)

 q = im[0:2047,0:767] - q1_dark
 valid = where(q1_mask gt 0 and q1_bpm lt 1)
 sky_vector = q[valid]
 mmm, sky_vector, skymod, skysigma, skyskew
 bad = where(q gt (skymod + sigma_dev*skysigma) and q1_mask gt 0, n_bad)
 for j = 0, n_bad -1 do begin
  value = skymod +  skysigma * RANDOMN(seed)
  q[bad[j]] = value
 endfor
 skycube[0:2047,0:767,ic] = q/median(q[valid])
 gains[0] = skymod
 print, skymod

 q = im[2048:4095,0:767] - q2_dark
 valid = where(q2_mask gt 0 and q2_bpm lt 1)
 sky_vector = q[valid]
 mmm, sky_vector, skymod, skysigma, skyskew
 bad = where(q gt (skymod + sigma_dev*skysigma) and q1_mask gt 0, n_bad)
 for j = 0, n_bad -1 do begin
  value = skymod +  skysigma * RANDOMN(seed)
  q[bad[j]] = value
 endfor
 skycube[2048:4095,0:767,ic] = q/median(q[valid])
 gains[1] = skymod
 print, skymod

 q = im[2048:4095,768:1535] - q3_dark
 valid = where(q3_mask gt 0 and q3_bpm lt 1)
 sky_vector = q[valid]
 mmm, sky_vector, skymod, skysigma, skyskew
 bad = where(q gt (skymod + sigma_dev*skysigma) and q1_mask gt 0, n_bad)
 for j = 0, n_bad -1 do begin
  value = skymod +  skysigma * RANDOMN(seed)
  q[bad[j]] = value
 endfor
 skycube[2048:4095,768:1535,ic] = q/median(q[valid])
 gains[2] = skymod
 print, skymod

 q = im[0:2047,768:1535] - q4_dark
 valid = where(q4_mask gt 0 and q4_bpm lt 1)
 sky_vector = q[valid]
 mmm, sky_vector, skymod, skysigma, skyskew
 bad = where(q gt (skymod + sigma_dev*skysigma) and q1_mask gt 0, n_bad)
 for j = 0, n_bad -1 do begin
  value = skymod +  skysigma * RANDOMN(seed)
  q[bad[j]] = value
 endfor
 skycube[0:2047,768:1535,ic] = q/median(q[valid])
 gains[3] = skymod
 print, skymod

 print, 'Processed '+ nam + '...'
endfor
writefits, tmp_path + 'skycube.fits', skycube


; Median sky
; ------------------------------
print, 'Now creating median sky and bad pixel map...'
sky = median(skycube,DIMENSION=3)
medsky = sky  ; save for tests below
writefits, out_path + 'sky.fits', sky*mask

; save gain values
; ----------------
gains = gains/mean(gains)
forprint, TEXTOUT= out_path + 'gains.txt', /NOCOMMENT, gains

; do comparison between simple median sky and masked sky
good = where(mask gt 0 and bpm lt 1)
sky = sky
sky[good] = sky[good]/flat[good]
sky = replace_pix(sky,xbad,ybad)
writefits, tmp_path + 'sky_nice.fits', sky * mask

smooth = filter_image(sky,FWHM=20,/ALL)
smooth_mask = filter_image(mask,FWHM=20,/ALL)
smooth[good] = smooth[good]/smooth_mask[good]
writefits, tmp_path + 'sky_smooth.fits', smooth * mask


print, 'sky.pro ended'

END
